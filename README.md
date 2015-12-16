# gatewayPrimerDesigner

Simple tool to attempt automation of gateway cloning primer design

##Required modules

The following perl modules are required by this script:

    Bio::Perl
    Bio::SeqIO
    Bio::DB::GenBank

These can be installed via cpan - see http://www.cpan.org/modules/INSTALL.html

In addition, perl version 5.10.0 or higher is required (run perl -v to see what version you have installed).

##Install and Run

Get the script:
    
    git clone https://github.com/gantzgraf/gatewayPrimerDesigner.git

Enter the newly downloaded directory:
 
    cd gatewayPrimerDesigner

Run the script to get usage information:
    
    perl gateway_designer.pl -h 

##Credit

Copyright © 2015  David A. Parry

Code for calculating primer TMs is lifted from perlprimer.pl (Copyright © 2003-2011, Owen Marshall)

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.
