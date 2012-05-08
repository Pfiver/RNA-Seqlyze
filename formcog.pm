#!/usr/bin/perl


 #***************************************************************
 # File:    formcog.pm
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

sub formcog {

    my ($pttfile1,$pttfile2,$circular1,$circular2,$allvsallfile,$genome2name,$bpgap,$genegap,$mapfile1,$mapfile2,$outfile,$eval)=@_;


    my @g1;
    my @ord1;
    my %id1;
    my $len1=read_genome_loc($pttfile1,\@g1,\@ord1,\%id1);
    
    my @g2;
    my @ord2;
    my %id2;
    my $len2=read_genome_loc($pttfile2,\@g2,\@ord2,\%id2);
    
    my %hit;
    process_hits($allvsallfile,$genome2name,\%id1,\%id2,$eval,\%hit);


    my @revord1;
    remap_genome(\@g1,$circular1,$bpgap,$genegap,$mapfile1,$len1,\@ord1,\@revord1,\%hit);
    

    my @revord2;
    remap_genome(\@g2,$circular2,$bpgap,$genegap,$mapfile2,$len2,\@ord2,\@revord2,\%hit);
    
    my %marked;
    compute_cogs($outfile,\%hit,\%marked,\%id1,\%id2,\@g1,\@g2,\@revord1,\@revord2);
}

sub remap_genome {
    my ($g,$circ,$bpgap,$genegap,$mapfile,$len,$ord,$revord,$hit)=@_;

    open(O,">$mapfile");
    print O "GeneId\tNewStart\tNewEnd\tStart\tEnd\n";

    my $n=scalar(@{$g});
    
    # find gap only if genome is circular
    if($circ) {
	
	my $beg=0;
	my $end=$n-1;

	my $genedist=1;

	while($beg<$n && !$$hit{$$g[$$ord[$beg]][0]}) { $beg++;$genedist++;}
	while($end>=0 && !$$hit{$$g[$$ord[$end]][0]}) { $end--;$genedist++;}


	if($beg<$end) {
	
	    my $bpdist=$$g[$beg][1]+$len-$$g[$end][2];

	    # remap only if gap not found at the genome end
	    if($bpdist<=$bpgap && $genedist<=$genegap) {
		while($beg+1<$n && gooddistloc($bpgap,$genegap,$n,$beg,\$end,$g,$ord,$hit)) {
		    $beg=$end;
		}
		if($beg+1<$n && $end<$n) { # split between $beg $beg+1
		    for(my $i=0;$i<=$beg;$i++) {
			my $j=shift @{$ord};
			print O $$g[$j][0],"\t",$$g[$j][1]+$len,"\t",$$g[$j][2]+$len,"\t",$$g[$j][1],"\t",$$g[$j][2],"\n";
			$$g[$j][1]+=$len;
			$$g[$j][2]+=$len;
			push(@{$ord},$j);
		    }
		}
	    }
	}
    }
    close(O);


    for(my $i=0;$i<$n;$i++) {
	$$revord[$$ord[$i]]=$i;
    }

}

sub gooddistloc {
    my ($bpgap,$genegap,$n,$beg,$end,$g,$ord,$hit)=@_;
    $$end=$beg+1;
    
    my $genedist=1;
    while($$end<$n && !$$hit{$$g[$$ord[$$end]][0]}) { $$end++;$genedist++;}
	
    if($$end<$n && $genedist<=$genegap && $$g[$$ord[$$end]][2]-$$g[$$ord[$$end]][1]<=$bpgap) { return(1);}
    
    return(0);
}

sub process_hits {
    my ($hitsfile,$name,$id1,$id2,$eval,$h)=@_;

    my $start=0;

    open(F,$hitsfile);
    while(<F>) {
	chomp;
	if(substr($_,0,1) eq ">") {
	    if($start) {last;}
	    if(substr($_,1) eq $name) { $start=1;}
	}
	else {
	    if($start) {
		my @a=split;
		my @b=split(/\|/,$a[0]);
		my $gid1=$b[1];
		my @b=split(/\|/,$a[1]);
		my $gid2=$b[1];

		if((!$eval || $a[-2]<=$eval) && exists($$id1{$gid1}) && exists($$id2{$gid2})) {
		    push(@{$$h{$gid1}},$gid2);
		    push(@{$$h{$gid2}},$gid1);
		}
	    }
	}
    }
    close(F);

}

sub compute_cogs {
    
    my ($outfile,$hit,$marked,$id1,$id2,$g1,$g2,$revord1,$revord2)=@_;

    my $cogno=1;

    my $fO;

    open($fO,">$outfile");

    print $fO "FAMILY_ID\tCHROMOSOME1\tCHROMOSOME2\n";

    my $id;
    foreach $id (keys %{$hit}) {
	if(!$$marked{$id}) {
	    my @cog1=();
	    my @cog2=();
	    my @list;
	    push(@list,$id);
	    if($$id1{$id}) { fill_cogs(1,\@cog1,\@cog2,$hit,$marked,@list);}
	    else {fill_cogs(0,\@cog1,\@cog2,$hit,$marked,@list);}

	    print_cogs($cogno,\@cog1,\@cog2,$fO,$id1,$id2,$g1,$g2,$revord1,$revord2);
	    $cogno++;
	}
    }
}

sub print_cogs {
    my ($cogno,$cog1,$cog2,$fh,$id1,$id2,$g1,$g2,$revord1,$revord2)=@_;

    print $fh "COG",$cogno,"\t";
    
    # print genes in genome 1
    for(my $i=0;$i<scalar(@{$cog1});$i++) {
	my $id=$$cog1[$i];

	if($i) { print $fh ":";}
	my $pos=$$id1{$id};

	print $fh $$g1[$pos][1],"-",$$g1[$pos][2],"-",$$revord1[$pos];
    }
    print $fh "\t";

    # print genes in genome 2
    for(my $i=0;$i<scalar(@{$cog2});$i++) {
	my $id=$$cog2[$i];

	if($i) { print $fh ":";}
	my $pos=$$id2{$id};

	print $fh $$g2[$pos][1],"-",$$g2[$pos][2],"-",$$revord2[$pos];
    }
    print $fh "\n";
}


sub fill_cogs {
    my ($genome,$cog1,$cog2,$hit,$marked,@list)=@_;

    while(@list) {
	my $id=pop(@list);
	if(!$$marked{$id}) {
	    $$marked{$id}=1;
	    if($genome) {
		push(@{$cog1},$id);
	    }
	    else { 
		push(@{$cog2},$id);
	    }
	    if(@{$$hit{$id}}) { fill_cogs(($genome+1)%2, $cog1, $cog2,$hit,$marked,@{$$hit{$id}});}
	}
    }
}

sub read_genome_loc {
    my ($pttfile,$g,$ord,$id)=@_;
    
    my $len;

    my $n=0;    
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

    my @unord;
    while(<F>) {
	my @rec;
	my @a=split;
	$rec[0]=$a[3]; # gi pid
	$$id{$a[3]}=$n;
	my @b=split(/\./,$a[0]);
	$rec[1]=$b[0];
	$rec[2]=$b[-1];
	$$id{$rec[0]}=$n;
	push(@{$$g[$n]},@rec);
	$unord[$n]=$n;
	$n++;
    }
    close(F);

    @{$ord} = sort bybeg @unord;
    
    sub bybeg {
	$$g[$a]->[1] <=> $$g[$b]->[1];
    }
    
    return($len);

}

