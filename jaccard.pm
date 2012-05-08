#!/usr/bin/perl



 #***************************************************************
 # File:    jaccard.pm
 #
 # Copyright (C) 2008  Mihaela Pertea (mpertea@umiacs.umd.edu)
 #
 # This program is free software; you can redistribute it and/or modify
 # it under the terms of the GNU General Public License as published by
 # the Free Software Foundation; either version 2 of the License, or
 # (at your option) any later version.
 #
 # This program is distributed in the hope that it will be useful,
 # but WITHOUT ANY WARRANTY; without even the implied warranty of
 # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 # GNU General Public License for more details.
 #
 # You should have received a copy of the GNU General Public License
 # along with this program; if not, write to the Free Software
 # Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 #
 #****************************************************************

use strict;

return 1;

sub jaccard_dist {
    
    my ($blastdir,$ngenomes, $namedir,$genes,$startdir,$color,$eval,$jdist,$maxjdist)=@_;

    my @hom;
    
    for(my $i=1;$i<=$ngenomes;$i++) {
	
	print STDERR "$i: ",$$namedir[$i],"\n";

	onegenome_compare($blastdir,$i,$$namedir[$i],$startdir,$color,$eval,\@hom);
    }


    for(my $i=1;$i<=$ngenomes;$i++) {
	$$maxjdist[$i-1]=0;
	for(my $j=1;$j<=$ngenomes;$j++) {
	    if($i!=$j) {
		$$jdist[$i-1][$j-1]=($$genes[$i]+$$genes[$j]-$hom[$i][$j]-$hom[$j][$i])/($$genes[$i]+$$genes[$j]);
		if($$jdist[$i-1][$j-1]>$$maxjdist[$i-1]) { $$maxjdist[$i-1]=$$jdist[$i-1][$j-1];}
	    }
	}
    }

}

sub how_many_genes {
    my ($filename)=@_;
    open(I,$filename);
    while(<I>) {
	chomp;
	my @a=split;
	if($a[-1] eq "proteins") {
	    close(I);
	    return($a[0]);
	}
    }
    close(I);
    print STDERR "No protein number found in the ptt file: ",$filename,"\n";
    return(0);
}

sub onegenome_compare {
    my ($blastdir,$n1,$name1,$startdir,$color,$eval,$hom)=@_;
    
    my $blastfile=$blastdir."/".$name1.".blast";

    my %seen;

    my $name2;
    my $n2;

    open(I,$blastfile);
    while(<I>) {
	chomp;
	my @a=split;
	if(substr($_,0,1) eq '>') {
	    $name2=substr($_,1);
	    $n2=$$color[$$startdir{$name2}];  # need to skip a genome if not in the ptt.list
	    %seen=();
	}
	else {
	    if($n2) {
		if($a[-2]<=$eval) {
		    my $gene=$a[0];
		    if(!$seen{$gene}) {
			$seen{$gene}=1;
			$$hom[$n1][$n2]++;
		    }
		}
	    }
	}
    }
    close(I);
}		    
