#!/usr/bin/perl
use strict;
use warnings; 

use Getopt::Long;
use Bio::Perl;
use Bio::SeqIO;
use Bio::DB::GenBank;
use List::Util qw (min);
use Data::Dumper;
use Term::ANSIColor;# qw(:constants);

my %codons = 
(
    AAA=>"K", AAC=>"N", AAG=>"K", AAT=>"N", 
    ACA=>"T", ACC=>"T", ACG=>"T", ACT=>"T", 
    AGA=>"R", AGC=>"S", AGG=>"R", AGT=>"S", 
    ATA=>"I", ATC=>"I", ATG=>"M", ATT=>"I", 
    CAA=>"Q", CAC=>"H", CAG=>"Q", CAT=>"H", 
    CCA=>"P", CCC=>"P", CCG=>"P", CCT=>"P", 
    CGA=>"R", CGC=>"R", CGG=>"R", CGT=>"R", 
    CTA=>"L", CTC=>"L", CTG=>"L", CTT=>"L", 
    GAA=>"E", GAC=>"D", GAG=>"E", GAT=>"D", 
    GCA=>"A", GCC=>"A", GCG=>"A", GCT=>"A", 
    GGA=>"G", GGC=>"G", GGG=>"G", GGT=>"G", 
    GTA=>"V", GTC=>"V", GTG=>"V", GTT=>"V", 
    TAA=>"*", TAC=>"Y", TAG=>"*", TAT=>"Y", 
    TCA=>"S", TCC=>"S", TCG=>"S", TCT=>"S", 
    TGA=>"*", TGC=>"C", TGG=>"W", TGT=>"C", 
    TTA=>"L", TTC=>"F", TTG=>"L", TTT=>"F"
);

my %opts = ();
GetOptions(
    \%opts,
    "c|cterm",  
    "colour|color",
    "d|number_translation",
    "g|gene=s",
    "h|?|help",
    "l|line_length=i",
    "max_tm=f",
    "min_tm=f",
    "n|nterm",
    "q|query=s",
    "v|verbose",
) or usage("Syntax error.\n");
usage() if $opts{h};
my $query = $opts{g} || ''; 
if ($opts{q}){
    $query .= "$opts{q}";
}elsif($opts{g}){
    $query .= " [GENE] AND Human [ORGN] and 0:10000 [SLEN]"; 
}
$opts{max_tm} ||= 90;
$opts{min_tm} ||= 45;
$opts{l} ||= 60;

my %primer_starts = 
(
    native_n => 'ggggacaagtttgtacaaaaaagcaggcttcgaaggagatagaaccatgg',
    fusion_n => 'ggggacaagtttgtacaaaaaagcaggcttc',
    native_c => 'ggggaccactttgtacaagaaagctgggtccta',
    fusion_c => 'ggggaccactttgtacaagaaagctgggtc',
);

if (not $query){
    usage("-g/--gene and/or  -q/--query options are required.\n");
}


my $query_obj = Bio::DB::Query::GenBank->new
( 
    -db => 'nucleotide', 
    -query => $query 
);


my $stream_obj = Bio::DB::GenBank->new->get_Stream_by_query($query_obj);

while (my $seq_object = $stream_obj->next_seq ) {
    #print Dumper $seq_object;
    processSeqObject($seq_object);
}

$stream_obj->DESTROY;
$query_obj->DESTROY;


###########################################################
sub processSeqObject{
    my ($seq_object) = @_;
    my %acc = ();
    for my $feat_object ( $seq_object->get_SeqFeatures ) {
        $acc{$seq_object->accession_number} = undef;
    }
    if (keys %acc > 1){
        die "More than one accession:\n" . join("\n", keys %acc) . "\nExiting\n";
    }
    my $accession = (keys %acc)[0];
    my @gene_features =
      grep { $_->primary_tag eq 'gene' } $seq_object->get_SeqFeatures;
    my @cds_features =
      grep { $_->primary_tag eq 'CDS' } $seq_object->get_SeqFeatures;
    my %cdna_seqs = 
      map { $_->get_tag_values('protein_id') => $_->spliced_seq->seq }
      @cds_features;
    my %protein_to_gene =
      map { $_->get_tag_values('protein_id') => $_->get_tag_values('gene') }
      @cds_features;
    my %gene_seqs = 
      map {$_->get_tag_values('gene'), $_->spliced_seq->seq } 
      @gene_features;
    if (keys %gene_seqs > 1){
        die "More than one gene sequence:\n" . 
        join("\n", keys %gene_seqs) . "\nExiting\n";
    }
    my $gene_sequence = $gene_seqs{(keys %gene_seqs)[0]};
    foreach my $cds (keys %cdna_seqs){
        my $identifier = "$accession/$cds";
        my $primer_pair = 0;
        my @hits = searchSequence
        (
            $gene_sequence, 
            $cdna_seqs{$cds}, 
        );
        if (not @hits){
            warn "Could not map CDS to gene sequence for '$cds'\n";
            next;
        }elsif (@hits > 1){
            warn "More than one mapping position for CDS in gene sequence for ".
                 "'$cds' - skipping\n";
            next;
        }
        my $l = length($cdna_seqs{$cds}); 
        foreach my $hit (@hits){
            my $start = $hit - length($cdna_seqs{$cds});
            my $end   = $hit;
            my @fseqs = ();
            my @rseqs = ();
            if ($opts{n}){ #n-terminal fusion 
                push @fseqs, substr($cdna_seqs{$cds}, 0, 25);
                push @fseqs, substr($cdna_seqs{$cds}, 3, 25);
            }else{#no n-terminal fusion - using Kozack seq of ATGG
                push @fseqs, substr($cdna_seqs{$cds}, 1, 25);
                push @fseqs, substr($cdna_seqs{$cds}, 4, 25);
            }
            if ($opts{c}){ #c-terminal fusion, avoid STOP
                push @rseqs, reverseComplement
                (
                    substr($cdna_seqs{$cds}, $l - 28, 25)
                );
            }else{#no c-terminal fusion - using endogenous STOP?
                push @rseqs,reverseComplement 
                (
                    substr($cdna_seqs{$cds}, $l - 25)
                );
            }
            my %f_prime = getPrimersAndTms(\@fseqs, $gene_sequence);
            if (not keys %f_prime){
                warn "No unique forward primers found within ". 
                "$opts{min_tm} < TM < $opts{max_tm} for $identifier\n";
                next;
            }
            my %r_prime = getPrimersAndTms(\@rseqs, $gene_sequence);
            if (not keys %r_prime){
                warn "No unique reverse primers found within ". 
                "$opts{min_tm} < TM < $opts{max_tm} for $identifier\n";
                next;
            }
            my @closest_pairs = getClosestTmPrimers(\%f_prime, \%r_prime);
            foreach my $pairTms (@closest_pairs){
                my @fseqs = @{$f_prime{$pairTms->[0]}};
                my @rseqs = @{$r_prime{$pairTms->[1]}};
                foreach my $f (@fseqs){
                    my $full_f = $opts{n} 
                      ? $primer_starts{fusion_n} 
                      : $primer_starts{native_n};
                    $full_f .= $f;
                    foreach my $r (@rseqs){
                        $primer_pair++;
                        my $full_r = $opts{c} 
                          ? $primer_starts{fusion_c} 
                          : $primer_starts{native_c};
                        $full_r .= $r;
                        my $header =  <<EOT

$identifier pair $primer_pair
Forward: $full_f
Reverse: $full_r
Targetting primers: $f/$r (TM $pairTms->[0]/$pairTms->[1])
EOT
;
                        if ($opts{colour}){
                            print color('bold');
                            print $header;
                            print color('reset');
                        }else{
                            print $header;
                        }
                        printAlignedCdsAndPrimers
                        (
                            $gene_sequence, 
                            $start,
                            $end,
                            $f,
                            $r,
                        );
                        printTranslation
                        (
                            $gene_sequence, 
                            $f,
                            $r,
                        );
                    }
                }
            }
        }
        
    }
    

}

###########################################################
sub printTranslation{
    my ($seq, $f, $r) = @_;
    my $full_f;
    my $full_r;
    my $translation_start = 0;
    if ($opts{n}){
        $full_f = "$primer_starts{fusion_n}$f";
        $translation_start = 4;
    }else{
        $full_f = "$primer_starts{native_n}$f";
        $translation_start = 46;
    }
    if ($opts{n}){
        $full_r = "$primer_starts{fusion_c}$r";
    }else{
        $full_r = "$primer_starts{native_c}$r";
    }
    $r = reverseComplement($r);
    $full_r = reverseComplement($full_r);
    my ($f_hit) = searchSequence($seq, $f);
#we've already checked that these are unique hits in getPrimersAndTms sub
    my ($r_hit) = searchSequence($seq, $r);
    my $clone_seq = $full_f . substr($seq, $f_hit, ($r_hit - length($r)) - $f_hit) . $full_r;
    my $peptide = translateOrf( substr($clone_seq, $translation_start)); 
    lineUpAndPrintTranslation($clone_seq, $peptide, $translation_start);     
}

###########################################################
sub translateOrf{
    my $dna = shift;
    my $protein = ''; 
    while ($dna =~ /(...)/g){
        my $codon = uc($1);
        my $aa = '?';
        if (exists $codons{$codon}){
            $aa = $codons{$codon};
        }
        $protein .= $aa;
        last if $aa eq '*';
    }
    return $protein; 
}

#########################################################
sub lineUpAndPrintTranslation{
    my ($dna, $protein, $start) = @_;
    $protein = join('', map { "-$_-" } split('', $protein) );
    $protein = '-' x $start . $protein;
    (my $raw_protein = $protein) =~ s/\-//g;
    my $protein_length = length($raw_protein); 
    my $dna_length = length($dna);
    my $header = <<EOT
PCR Product Length: $dna_length
Translation Length: $protein_length

EOT
;

    if ($opts{colour}){
        print color('bold');
        print $header;       
        print color('reset');
    }else{
        print $header;       
    }
    #my $primer_line = " " x length($coding_upper); 
    my $num_length = length($dna_length); 
    for (my $i = 0; $i < $dna_length; $i += $opts{l}){
        my $l = $opts{l};
        if ($i + $opts{l} > $dna_length){
            $l = $dna_length - $i ;
        }
        my $pl = $opts{l};
        if ($i + $opts{l} > length($protein)){
            $pl = length($protein) - $i ;
        }
        my $d = substr($dna, $i, $l);
        my $p = '';
        if ($i < length($protein)){
            $p = substr($protein, $i, $pl);
        }
        if ($opts{d}){
            my $n = $i + 1; #dna position
            $d = sprintf("%${num_length}d: %s", $n, $d);
            my $cn = $i - $start; #0-based coding DNA position of first letter of line
            my $pn = 1;#protein position... 
            if ($i > 0){
                $pn = int($cn/3) + 1 ;
                $pn++ if ($cn % 3 == 2); #if first nucleotide of line is the last letter of 
                                     #codon we represent the next protein 'letter'
            }
            $pn = $pn <= $protein_length ? $pn : $protein_length;
            $p = sprintf("%${num_length}d: %s", $pn, $p);
            
        }
        print "$d\n$p\n\n";
    } 
}
###########################################################
sub printAlignedCdsAndPrimers{
    my (
        $seq, 
        $cds_start, 
        $cds_end, 
        $fprime, 
        $rprime
       ) = @_;
    $rprime = reverseComplement($rprime);
    my $coding_upper = lc (substr($seq, 0, $cds_start));
    $coding_upper .= uc(substr($seq, $cds_start, $cds_end - $cds_start));
    $coding_upper .= lc(substr($seq, $cds_end - $cds_start));
    my ($f_hit) = searchSequence($seq, $fprime);
#we've already checked that these are unique hits in getPrimersAndTms sub
    my ($r_hit) = searchSequence($seq, $rprime);
    my $f_start = $f_hit - length($fprime);
    my $r_start = $r_hit - length($rprime);
    my $target_length = $r_hit - $f_start;
    if ($opts{colour}){
        print color('bold');
        print "Target Length: $target_length\n";
        print color('reset');
    }else{
        print "Target Length: $target_length\n";
    }
    #my $primer_line = " " x length($coding_upper); 
    my $arrow_line = " " x length($coding_upper); 
    #substr($primer_line, $f_start, length($fprime), $fprime);
    #substr($primer_line, $r_start, length($rprime), $rprime);
    substr($arrow_line, $f_start , length($fprime), ">" x length($fprime));
    substr($arrow_line, $r_start , length($rprime), "<" x length($rprime));
    for (my $i = 0; $i < length($coding_upper); $i += $opts{l}){
        my $l = $opts{l};
        if ($i + $opts{l} > length($coding_upper)){
            $l = length($coding_upper) - $i ;
        }
        my $s = substr($coding_upper, $i, $l) ;
        if ($opts{colour}){
            my @color_coords = ();
            push @color_coords, getPrimerColorCoords($i, $l, $r_start, $r_start + length($rprime));
            push @color_coords, getPrimerColorCoords($i, $l, $f_start, $f_start + length($fprime));
            for my $ar (@color_coords) {
                my ($start, $len) = @$ar;
                $s = substr($s, 0, $start)                  # first part
                   . colored(substr($s,$start,$len), 'red' )  # colored part
                   . substr($s,$start+$len);                  # final part
            }
        }
        print "$s\n"; 
            #my $pr = substr($primer_line , $i, $l) . "\n"; 
        my $ar =  substr($arrow_line , $i, $l) . "\n"; 
        if ($ar =~ /\S/){
            #print $pr;
            print $ar;
        }else{
            print "\n";
        }
    }
    print "\n";
        
}

###########################################################
sub getPrimerColorCoords{
    my ($offset, $length, $start, $end) = @_;
    return if $start > $offset + $length;
    return if $end < $offset; 
    my $s = $start - $offset; 
    my $e = $end - $offset;
    $s = $s > 0 ? $s : 0;
    $e = $e < $length ? $e : $length ;
    return [$s,  $e - $s];#return anon array of index and length
}

###########################################################
sub reverseComplement{
    my $seq = shift;
    $seq = reverse($seq);
    $seq =~ tr/acgtACGT/tgcaTGCA/;
    return $seq;
}

###########################################################
sub getClosestTmPrimers{
    #takes two hashes of primer TMs to ref to array of primer seqs
    #returns array of array refs to pairs primer TMs with closest values
    my ($fs, $rs) = @_;
    my @ftms = sort {$a <=> $b} keys %$fs;
    my @rtms = sort {$a <=> $b} keys %$rs;
    my %diff_pairs = ();
    my $prev_diff; 
FOR: for my $f (@ftms){
        my $r_diff;
REV:    for my $r (@rtms){
            my $diff = $f - $r;
            if (not defined $prev_diff or abs($diff) <= abs($prev_diff)){
                push @{$diff_pairs{$diff}}, [$f, $r];
                $prev_diff = $diff;
            }else{
                if (not defined $r_diff or abs ($diff) <= abs($r_diff)){
                    $r_diff = $diff;
                }else{#this diff is greater than our last in @rtms
                    last REV;
                }
            }
        } 
    } 
    my $min = min( map { abs($_) } keys (%diff_pairs) ) ;
    my @closest = ();
    push @closest, @{$diff_pairs{$min}} if exists $diff_pairs{$min};
    push @closest, @{$diff_pairs{$min*-1}} if exists $diff_pairs{$min*-1};
    return @closest;
}
###########################################################
sub getPrimersAndTms{
    #returns a hash of primer TMs to ref to array of primer seqs
    my $seqs = shift;
    my $dna = shift;
    my %tm_to_primers = (); 
    foreach my $seq (@$seqs){
        for my $length (18..25){
            for (my $i = $length; $i  <= length($seq); $i++){
                my $p = substr($seq, 0, $i);
                my @hits = searchSequence($dna, $p); 
                next if @hits > 1;#skip non-unique mapping primers
                my $tm = calculateTm($p);
                next if $tm > $opts{max_tm};
                next if $tm < $opts{min_tm};
                push @{$tm_to_primers{$tm}}, $p;
            }
        }
    }
    return %tm_to_primers;
}
###########################################################
sub calculateTm{
    my $dna = shift;
    $dna = uc($dna);
	my $at = 0;
	my $cg = 0;
	while ($dna =~ /(.)/g){
		if ($1 =~ /[AT]/){
			$at++;
		}elsif ($1 =~ /[CG]/){
			$cg++;
		}else{
			die "illegal character (\"$1\") in DNA sequence\n";
		}
	}
	return (2 * $at) + (4 * $cg);
}

###########################################################
sub searchSequence {
    my ( $dna, $search_sequence ) = @_;
    my @found = ();
    while ( $dna =~ /$search_sequence/ig ) {
        push( @found, pos($dna) );
    }
    return @found;
}

###########################################################

sub usage {
    my $msg = shift;
    if ($msg) {
        print STDERR "\n$msg";
    }
    
    print STDERR <<EOT

Options: 
    -g,--gene  <gene query to search genbank with>
        By default this searches for human sequences of less than 10,000 bp (i.e. gene [GENE] AND Human [ORGN] and 0:10000 [SLEN]). However, you may change the parameters used using the -q/--query option.
    -q,--query <Detailed query/parameters>
        Use instead of -g query if you wish to enter a detailed query or use in conjunction with -g/--gene option to specify custom parameters for the search. If you know the accession of your sequence you can enter it here for a quicker and more specific search.
    -c,--cterm
        Use this flag to design for C-Terminal fusion proteins.
    -n,--nterm
        Use this flag to design for N-Terminal fusion proteins.
    -d,--number_translation
        Use this flag to number the translated clone output.
    -l,--line_length
        Specify the length of lines for your sequence output. Default = 60.
    --max_tm
        Maximum TM for primers. Default = 90.
    --min_tm
        Minimum TM for primers. Default = 45.
    --colour,--color
        Use this flag to colour your output.
    -h,--help
        Show this message and exit.
EOT
;
    exit 1 if $msg;
    exit;
}
