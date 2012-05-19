#!/usr/bin/perl

#***************************************************************
 # File:    confidence.pm
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

sub print_operon {
    my ($fh,$operon,$pair,$conf)=@_;

    my $ncons=scalar(@{$operon})/2;
    
    printf $fh "%s %d %d",$pair,$conf,$ncons;
    for(my $j=0;$j<$ncons;$j++) {
	print $fh " ",$$operon[2*$j+1];
    }
    print $fh "\n";

}

sub compute_conf {  
    my ($genomen,$operon,$dist,$problen,$jdist,$maxjdist)=@_;
    
    my $n=scalar(@{$operon})/2;

    my $conf=0;

    for(my $j=0;$j<$n;$j++) {

#	print STDERR "Genome of consv pair is:",$$operon[2*$j],"\n";

	my $stat=$$dist[$genomen][$$operon[2*$j]][2];  # statistics are oSpairs = S pairs with orientation conserved
	my $bg=$$dist[$genomen][$$operon[2*$j]][0]-$stat+$$dist[$genomen][$$operon[2*$j]][1];    # background = Spairs-oSpairs+Dpairs

	if(!$stat) { print STDERR "Error - no S(d) conserved pairs for operon in $genomen with pair ",$$operon[2*$j+1],"!!!\n"; }
	else {
	    my $locconf=$$jdist[$genomen][$$operon[2*$j]]*(1-$problen*$bg/$stat)*100/$maxjdist;
	    if($locconf>$conf) { $conf=$locconf;}
	}

    }

    return ($conf);
}

sub read_pairs {
    my ($filename,$operon,$n,$consv,$maxdist,$g,$id,$len,$ind)=@_;

    open(F,$filename);
    while(<F>) {
	chomp;
	my @a=split;

	if($a[2] == $consv && gooddist($a[0],$a[1],$maxdist,$g,$id,$len,$ind)) {
	    my $pair=$a[0]." ".$a[1];
	    my @rec;
	    $rec[0]=$n;
	    $rec[1]=$a[3]." ".$a[4];
	    push(@{$$operon{$pair}},@rec);
	}
    }
    close(F);
}

sub gooddist {
    my ($id1,$id2,$maxdist,$g,$id,$len,$ind)=@_;
    
    my $n1=$$id{$id1};
    my $n2=$$id{$id1};

    if($n1>$n2) { 
	my $n=$n2;
	$n2=$n1;
	$n1=$n;
    }

    my $i=0;

    my $nind=scalar(@{$ind});

    while($i<$nind-1 && $n1>=$$ind[$i+1]) { $i++;}
    
    my $dist=$$g[$n2][1]-$$g[$n1][2];
    if($dist<$maxdist) { return(1);}

    $dist=$$g[$n1][1]+$$len[$i]-$$g[$n2][2];
    if($dist<$maxdist) { return(1);}

    return(0);
}

	  
sub read_conf {
    my ($filename)=@_;

    my ($conf,$dconf,$oconf,$odconf);

    open(I,$filename);
    while(<I>) {
	chomp;
	
	my @a=split;

	if($a[0] eq "Confidence") { 
	    $conf=$a[-1];
	}
	elsif($a[0] eq "Confidence(d)") {
	    $dconf=$a[-1];
	}
	elsif($a[0] eq "o-Confidence") { 
	    $oconf=$a[-1];
	}
	elsif($a[0] eq "o-Confidence(d)") {
	    $odconf=$a[-1];
	}
    }
    close(I);

    return($conf,$dconf,$oconf,$odconf);
}

sub read_genome {
    my ($pttfn,$pttfile,$g,$id,$n,$maxdist,$binlen,$poslen,$neglen)=@_;
    
    my $len;

    open(F,$pttfile);
    
    for(my $i=0;$i<3;$i++) {
	my $line;

	do {
	    $line=<F>;
	    chomp($line);
	} while(!$line);
	
	if(!$i) {
	    my @a=split(/\./,$line);
	    $len=$a[-1];
	}
    }

    my $lastsign="";
    my $lastend=0;

    my $maxbin=int($maxdist/$binlen);

    my $possum=0; # use pseudocounts of 1
    my $negsum=0;

    while(<F>) {
	my @rec;
	my @a=split;
	$rec[0]=$a[3]; # gi pid
	$$id{$a[3]}=$$n;
	my @b=split(/\./,$a[0]);
	$rec[1]=$b[0];
	$rec[2]=$b[-1];

	my $genedist=$rec[1]-$lastend;
	    
	if($genedist>0) { # just check the positive genes, overalping genes are estimated from non-overlaping ones due to other mechanisms for overlaping genes
	    my $bin=int($genedist/$binlen);
	    if($bin>$maxbin) { $bin=$maxbin;}
	    if($lastsign eq $a[1]) { $$poslen[$pttfn][$bin]++; $possum++;} # S pair
	    elsif($lastsign ne "") {$$neglen[$pttfn][$bin]++; $negsum++;} # D pair
	}

	$lastsign=$a[1];
	$lastend=$rec[2];

	$$id{$rec[0]}=$$n;
	push(@{$$g[$$n]},@rec);
	$$n+=1;
    }
    close(F);

    for(my $i=0;$i<=$maxbin;$i++) {
	$$poslen[$pttfn][$i]++; # use pseudocounts of 1
	$$neglen[$pttfn][$i]++;
	$$poslen[$pttfn][$i]/=$possum+$maxbin;
	$$neglen[$pttfn][$i]/=$negsum+$maxbin;
    }

    return($len);

}

sub intergenic_distance_prob {
    my ($operon,$g,$id,$binlen,$maxdist,$len,$ind,$poslen,$neglen)=@_;
    
    my $maxbin=int($maxdist/$binlen);

    my($id1,$id2)=split(/\s+/,$operon);

    my $n1=$$id{$id1};
    my $n2=$$id{$id1};

    if($n1>$n2) { 
	my $n=$n2;
	$n2=$n1;
	$n1=$n;
    }

    my $i=0;

    my $nind=scalar(@{$ind});

    while($i<$nind-1 && $n1>=$$ind[$i+1]) { $i++;}
    
    my $dist=$$g[$n2][1]-$$g[$n1][2];
    my $fardist=$$g[$n1][1]+$$len[$i]-$$g[$n2][2];
    if($fardist<$dist) { $dist=$fardist;}

    my $bin=0;
    if($dist>0) {
	$bin=int($dist/$binlen);
	if($bin>$maxbin) { $bin=$maxbin;}
    }

    my $prob=$$neglen[$i][$bin]/$$poslen[$i][$bin];

    return($prob);
}
