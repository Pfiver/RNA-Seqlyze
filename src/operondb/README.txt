OperonDB -- Dec 2008

Copyright (C) 2008 Mihaela Pertea (mpertea@umiacs.umd.edu)
---------------------------------------------------------------------
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or (at
your option) any later version.
 
This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
USA
---------------------------------------------------------------------

This README file explains how to use the script createoperondb.pl, and 
what preliminary steps are needed. 

1. First compile the HomologyTeam software by going to the homologyteams 
directory and simply typing 'make'. This should result in an executable
named 'homologyteams'. The software resulted is a slightly 
modified version of the software implemented in:

Identifying Conserved Gene Clusters in the Presence of Orthologous Groups
Xin He and Michael H. Goldwasser
Proceedings of the Eighth Annual International Conferences on Research in 
Computational Molecular Biology (RECOMB), San Diego, California, Mar. 2004, 
pp. 272-280.

For more information about the initial version of this software please
check  http://euler.slu.edu/~goldwasser/homologyteams/

2. Create a database of all prokaryotic genomes. We did this by downloading
all prokaryotic genomes available from Genbank (www.ncbi.nih.gov). 
 A unique identifier name should be assigned to each genome, and this name 
should be used to create a directory with the same name for each genome 
containing all data of interest for that genome. An example of such a 
directory for the Acinetobacter sp. ADP1 is provided in 
      examples/genomes/Acinetobacter_sp_ADP1
When downloading the prokaryotic genomes from Genbank we were particularly 
interested in the protein sequences of the genes in each genome (the *.faa
files), and the files containing the descriptions of all protein coding
genes (the *.ptt) files.
 Create a file with the location of all *.ptt files. An example of such a file is
given in 
      examples/ptt.list

3. You will need to run blast searches of all gene coding sequences in a 
prokaryotic genome against all gene coding sequences in all other 
prokaryotic genomes. 
 The blast searches corresponding to a prokaryotic genome should be put in a 
file named with the unique name assigned to that genome followed by the 
extension ".blast".  All *.blast files should be put in the same directory.
 When creating the searches for OperonDB we used the blastall command with 
the following parameters:
 blastall -p blastp -e 1e-05 -m 8 -v 0 -F f 
on all *.faa files.
 The blast file corresponding to a given genome should be organized like this: 
all searches between that genome and another one should be preceded by the 
symbol ">" followed by the name of the other genome. An example file containing
all searches of all genes in the genome named Acinetobacter_sp_ADP1 against all
other genomes in the database is given in 
      examples/blastdir/Acinetobacter_sp_ADP1.blast

4. Now you are ready to run the createoperondb.pl script by using the following 
command: 
     createoperondb.pl <list_of_ptt_files> <directory_of_all_blast_searches> <homologyteam_program>
E.g. using the above examples the command above would look like this:
     createoperondb.pl examples/ptt.list examples/blastdir homologyteams/homologyteams	
Of course this will not work on our limited examples since the database is not 
complete, but hopefully it will work on the database you created. If not, please send 
your problems or comments to:
http://www.cbcb.umd.edu/cgi-bin/mailer_scripts/contact_us.pl?receiver=mpertea&refer_text=Operon
