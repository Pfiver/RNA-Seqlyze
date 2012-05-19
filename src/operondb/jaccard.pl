#!/usr/bin/perl


use strict;
use Getopt::Long;

use FindBin;
use lib "$FindBin::Bin";

use computedist;
use jaccard;

my ($pttlist,$blastdir)=@ARGV;

my $eval=1e-15;

GetOptions( "e:s"    => \$eval
	    );


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

for(my $i=1;$i<=$ngenomes;$i++) {
  for(my $j=1;$j<=$ngenomes;$j++) {
    if($i!=$j) {
      print "$i $j ",$namedir[$i]," ",$namedir[$j]," ",$jdist[$i-1][$j-1]," ",$maxjdist[$i-1],"\n";
    }
  }
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
