#!/usr/bin/perl


 #***************************************************************
 # File:    createoperondb.pl
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
 # this script takes as input blast searches of all genomes against 
 # all genomes and predicts all possible operons in all genomes in 
 # the database
 #****************************************************************

use strict;
use Getopt::Long;

use FindBin;
use lib "$FindBin::Bin";

use computedist;
use formcog;
use confidence;
use jaccard;

my ($pttlist,$blastdir,$homologyprogram)=@ARGV;

my $usage = q/
  Creates operon predictions for all genomes described by a set of *.ptt 
  files and blast searches of all genomes against all others in the database

  Usage:
    createoperondb.pl <list_of_ptt_files> <directory_of_all_blast_searches> <homologyteam_program> [optional_parameters]

  [optional_parameters]

  -b bpdist      gives bp distance allowed between two genes in a prediction (default 1000)
  -d genedist    allowed gene distance (default 2)
  -l binlen      bin length for estimating the length distributions (default 200)
  -c consv       conservation type allowed between genes:
                 * 1 if only orientation is conserved
                 * 2 if both orientation and order are conserved (default)
  -e eval        maximum E-val value for homology
  -h             help

/;

### read options

my $bpdist=1000;
my $genedist=2;
my $conservationtype=2;
my $eval=1e-15;
my $help=0;
my $binlen=200;

GetOptions( "h"      => \$help,
	    "b:i"    => \$bpdist,
	    "d:i"    => \$genedist,
	    "c:i"    => \$conservationtype,
	    "l:i"    => \$binlen,
	    "e:s"    => \$eval
	    );


die $usage if($help);

my $ngenomes=0;
my @name;
my %startdir;
my @namedir;
my @color;
my @ptt;
my @file;
my @circ;
my @genes;

### read all ptt files:

my $n=1;   # ---> counts all ptt files
my $lastdirname="";

open(F,$pttlist);
while(<F>) {
    chomp;
    my @a=split(/\//);  # $a[-1] gives the ptt file name; $a[-2] gives the genome dir name
    $name[$n]=$a[-2];
    if($name[$n] ne $lastdirname) { 
	$ngenomes++;
	$startdir{$name[$n]}=$n;
	$namedir[$ngenomes]=$name[$n];
    }
    $color[$n]=$ngenomes;     # colors all ptt files
    my @b=split(/\./,$a[-1]);
    $ptt[$n]=$b[0];
    $file[$n]=$_;
    $circ[$n]=read_circ($file[$n]);  # circularity might be better read from the genbank file (*.gbk) directly
    $lastdirname=$name[$n];

    $genes[$ngenomes]+=how_many_genes($file[$n]);

    $n++;
}
close(F);

### compute the Jaccard distance

my @jdist;
my @maxjdist;

print STDERR "Compute Jaccard distance for $ngenomes genomes:\n"; 

jaccard_dist($blastdir,$ngenomes, \@namedir,\@genes,\%startdir,\@color,$eval,\@jdist,\@maxjdist);

### search all genomes against all others

my $lastdref=0;
my $itg=1;

my @dist;
my %operon;
my @ind;
my @len;
my @g;
my %id;
my @poslen;
my @neglen;

print STDERR "Predict operons:\n";

while($itg<=$ngenomes*$ngenomes) {
    
    my $dref=1+int(($itg-1)/$ngenomes);  # genome of reference

    if($dref!=$lastdref) { # last genome of reference has been processed

	# compute operons
	if($lastdref) { # only if not first genome

	    my $operonoutfile=$namedir[$lastdref].".operons";
	
	    open(O,">$operonoutfile");
	    foreach my $o (keys %operon) {

		my $conf;
		my $problen=intergenic_distance_prob($o,\@g,\%id,$binlen,$bpdist,\@len,\@ind,\@poslen,\@neglen);


		my $conf=compute_conf($lastdref-1,\@{$operon{$o}},\@dist,$problen,\@jdist,$maxjdist[$lastdref-1]); # n gives the number of genomes where the operon is conserved

		print_operon(*O,\@{$operon{$o}},$o,$conf); # print the most conservative confidence (minimum of all 4)
		
	    }
	    close(O);

	    print STDERR "Done genome $lastdref with name ",$namedir[$lastdref],"\n";

	    # clean memory
	    @dist=();
	    %operon=();
	    @ind=();
	    @len=();
	    @g=();
	    %id=();
	    @poslen=();
	    @neglen=();
	}

	# process new genome
	my $nref=$startdir{$namedir[$dref]};
	my $colorref=$color[$nref];
	my $i=$nref;
	my $n=0;
	while($color[$i]==$colorref) {
	    $ind[$i-$nref]=$n;
	    $len[$i-$nref]=read_genome($i-$nref,$file[$i],\@g,\%id,\$n,$bpdist,$binlen,\@poslen,\@neglen);
	    $i++;
	}

    }
    else {

	my $it = 1+ ($itg-1) % $ngenomes;    # genome against which pairs are computed


	my $nref=$startdir{$namedir[$dref]};
	my $genomedir=$name[$nref];   
	my $blastfile=$blastdir."/".$genomedir.".blast";
	die "No file $blastfile!" if(!(-e $blastfile));


	# here consider all pairs of files from genomedir and namedir[$it]
	my $colorref=$color[$nref];
	if($it != $colorref) {
	    my $colorque=$it;
	    my $nque=$startdir{$namedir[$colorque]};
	    
	    my $i=$nref;
	    my $colorref=$color[$nref];

	    my ($ngenes,$ndir,$sdpairs,$ddpairs,$spairs,$dpairs,$csdpairs,$cddpairs,$cspairs,$cdpairs,$ocsdpairs,$ocddpairs,$ocspairs,$ocdpairs);

	    my %seen;
	    my %pairs;


	    while($color[$i]==$colorref) {
	    

		my $added=0;
		my $j=$nque;

		while($color[$j]==$colorque) {
	    

		    die "No ptt file:".$file[$i]." ".$file[$j] if(!(-e $file[$i]) || !(-e $file[$j]));
		    
		    formcog($file[$i],$file[$j],$circ[$i],$circ[$j],$blastfile,$name[$j],$bpdist,$genedist,$ptt[$i]."-".$ptt[$j].".map",$ptt[$j]."-".$ptt[$i].".map",$ptt[$i]."-".$ptt[$j].".cog",$eval);


		    my $cogfile=$ptt[$i]."-".$ptt[$j].".cog";
		    open(F,$cogfile);
		    <F>;
		    if(<F>) {

			my $command2="$homologyprogram -b 2 -d $bpdist -n $genedist -W ".$ptt[$i]."-".$ptt[$j].".out -O witness ".$ptt[$i]."-".$ptt[$j].".cog";
			system($command2);
		    

			my ($tngenes,$tndir,$tsdpairs,$tddpairs,$tspairs,$tdpairs,$tcsdpairs,$tcddpairs,$tcspairs,$tcdpairs,$tocsdpairs,$tocddpairs,$tocspairs,$tocdpairs)=computedist($ptt[$i]."-".$ptt[$j].".out",$file[$i],$file[$j],$ptt[$i]."-".$ptt[$j].".map",$ptt[$j]."-".$ptt[$i].".map",$circ[$i],$bpdist,$genedist,\%seen,\%pairs);
		    
		
			if(!$added) {
			    $ngenes+=$tngenes;
			    $ndir=$tndir;
			    $sdpairs+=$tsdpairs;
			    $ddpairs+=$tddpairs;
			    $spairs+=$tspairs;
			    $dpairs+=$tdpairs;
			    $added=1;
			}
			
			$csdpairs+=$tcsdpairs;
			$cddpairs+=$tcddpairs;
			$cspairs+=$tcspairs;
			$cdpairs+=$tcdpairs;
			$ocsdpairs+=$tocsdpairs;
			$ocddpairs+=$tocddpairs;
			$ocspairs+=$tocspairs;
			$ocdpairs+=$tocdpairs;
			
			my $command4="rm ".$ptt[$i]."-".$ptt[$j].".out";
			system($command4);

		    }
		    close(F);
		    my $command1="rm ".$ptt[$i]."-".$ptt[$j].".map ".$ptt[$j]."-".$ptt[$i].".map ".$ptt[$i]."-".$ptt[$j].".cog";
		    system($command1);
		    
		    
		    $j++;
		}

		$i++;
	    }


	    # option 1  --- before
	    #@{$dist[$dref-1][$it-1]}=($cspairs,$cdpairs,$ocspairs);

	    # option 2  --- plus non-adjacent genes 
    	    @{$dist[$dref-1][$it-1]}=($csdpairs,$cddpairs,$ocsdpairs);

	    # read potential operons

	    foreach my $k (keys %seen) {
		my ($g1,$sign1,$g2,$sign2)=split(/\./,$k);


		if($sign1==$sign2 && $seen{$k}==$conservationtype && gooddist($g1,$g2,$bpdist,\@g,\%id,\@len,\@ind)) {
		    my $pair=$g1." ".$g2;
		    my @rec;
		    $rec[0]=$it-1;
		    $rec[1]=$pairs{$k};
		    push(@{$operon{$pair}},@rec);
		}
	    }

	}
    }

    $itg++;
    $lastdref=$dref;
}

if($lastdref) { # only if not first genome

    my $operonoutfile=$namedir[$lastdref].".operons";
	
    open(O,">$operonoutfile");
    foreach my $o (keys %operon) {
	
	my $conf;
	my $problen=intergenic_distance_prob($o,\@g,\%id,$binlen,$bpdist,\@len,\@ind,\@poslen,\@neglen);


	my $conf=compute_conf($lastdref-1,\@{$operon{$o}},\@dist,$problen,\@jdist,$maxjdist[$lastdref-1]); # n gives the number of genomes where the operon is conserved

	print_operon(*O,\@{$operon{$o}},$o,$conf); # print the most conservative confidence (minimum of all 4)
	
    }
    close(O);

    print STDERR "Done genome $lastdref with name ",$namedir[$lastdref],"\n";
}


sub multiplier {
    my ($spairs,$dpairs,$cspairs,$cdpairs,$ngenes,$ndir,$genegap)=@_;
    
    my $multiplier=0;
    my $fraction=0;
    my $prob=0;
    my $confidence=0;

    if($spairs) {
	$multiplier-=$genegap*$ngenes;
	$multiplier+=2*$genegap*$ndir+$spairs;
	$multiplier/=$spairs;
    }
    else {
	print STDERR "No S pairs!\n";
    }

    if($dpairs && $cspairs) { 
	$fraction=$cdpairs*$spairs/($cspairs*$dpairs);
    }
    $prob=$fraction*$multiplier/2;
    $prob=1-$prob;
    $confidence=1-$fraction;
    return($multiplier,$fraction,$prob,$confidence);
}

sub read_circ {
    my ($filename)=@_;

    open(I,$filename);

    my $line="";
    do {
	$line=<I>;
	chomp($line);
    } while(!$line);
    
    close(I);
    
    if(index($line,"linear")>-1) {
	return(0);
    }
    else {
	return(1);
    }
}


