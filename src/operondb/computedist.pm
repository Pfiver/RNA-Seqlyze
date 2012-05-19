#!/usr/bin/perl


 #***************************************************************
 # File:    computedist.pm
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

sub computedist {
    my ($cogoutfile,$pttfile1,$pttfile2,$mapfile1,$mapfile2,$circular1,$bpgap,$genegap,$seen,$pairs)=@_;

    my %map1;
    my %map2;

    readmap($mapfile1,\%map1);
    readmap($mapfile2,\%map2);

    my @g1;
    my %id1;
    my $shift=read_genome_cd($pttfile1,\@g1,\%map1,\%id1);

    my ($sdpairs,$ddpairs,$spairs,$dpairs,$ngenes,$ndir)=compute_mult(\@g1,$circular1,$bpgap,$genegap,$shift);

    my @g2;
    my %id2;
    read_genome_cd($pttfile2,\@g2,\%map2,\%id2);

    my ($csdpairs,$cddpairs,$cspairs,$cdpairs,$ocsdpairs,$ocddpairs,$ocspairs,$ocdpairs)=compute_conserved($cogoutfile,\@g1,\@g2,\%id1,\%id2,$seen,$pairs,$genegap);

    return($ngenes,$ndir,$sdpairs,$ddpairs,$spairs,$dpairs,$csdpairs,$cddpairs,$cspairs,$cdpairs,$ocsdpairs,$ocddpairs,$ocspairs,$ocdpairs);
}


sub read_cogs {
    my ($line,$cog,$id)=@_;

    chomp($line);
    my $n=0;

    my @b=split(/\,/,$line);
    for(my $i=0;$i<=$#b;$i++) {
	chop($b[$i]);
	my @c=split(/:/,$b[$i]);
	
	my @rec;
	$rec[0]=$$id{$c[1]};
	my $pos=index($c[0],"(");
	$rec[1]=substr($c[0],$pos+1);

	push(@{$$cog[$n]},@rec);
	$n++;
    }
    
    return($n);
}

sub read_map_cogs {
    my ($line,$cog1,$cog2,$g,$id)=@_;
    chomp($line);

    my @b=split(/\,/,$line);
    for(my $i=0;$i<=$#b;$i++) {
	chop($b[$i]);
	my @c=split(/:/,$b[$i]);

	my $pos=index($c[0],"(");
	my $cogname=substr($c[0],$pos+1);
	my $n=$$id{$c[1]};
	push(@{$$cog2{$cogname}},$n);

	if(exists $$cog1{$cogname}) {
	    if($$cog1{$cogname}) {
		if($$g[$n][1]!=$$cog1{$cogname}) {
		    $$cog1{$cogname}=0;
		}
	    }
	}
	else {
	    $$cog1{$cogname}=$$g[$$id{$c[1]}][1];
	}
    }
}


sub compute_conserved {
    my ($cogfile,$g1,$g2,$id1,$id2,$seen,$pairs,$genegap)=@_;

    my $csdpairs=0;  ### conserved at gap distance
    my $cddpairs=0;
    my $cspairs=0;   ### conserved and adjacent
    my $cdpairs=0;
    my $ocsdpairs=0;  ### conserved at gap distance and keeping order -A-> -B-> : -A'-> -B'-> or <-B'- <-A'-
    my $ocddpairs=0;  # -A-> <-B- : -A'-> <-B'- or -B'-> <-A'- ; and <-A- -B-> : <-A'- -B'-> or <-B'- -A'-> 
    my $ocspairs=0;   ### conserved and adjacent and keeping order
    my $ocdpairs=0;
    my $ncog1;
    my @cog1;

    open(F,$cogfile);
    while(<F>) {
	chomp;
	my %cog2;
	my %cog3;
	if($_ eq "Witness in the 1st chromosome:") {
	    my $line=<F>;
	    @cog1=();
	    $ncog1=read_cogs($line,\@cog1,$id1);
	}
	elsif($_ eq "Witness in the 2nd chromosome:") {
	    my $line=<F>;
	    %cog2=();
	    %cog3=();
	    read_map_cogs($line,\%cog2,\%cog3,$g2,$id2);
	    for(my $i=0;$i<$ncog1-1;$i++) {

		my $n1=$cog1[$i][0];
		my $j=$i+1;
		my $n2=$cog1[$j][0];
		my $s1=$cog2{$cog1[$i][1]};

		my $pair=$$g1[$n1][3].".".$$g1[$n1][1].".".$$g1[$n2][3].".".$$g1[$n2][1];

		# $seen{$pair} = 1 if the pair was encountered before *and* there was orientation conservation observed, but not order conservation
		# $seen{$pair} = 2 if the pair was encountered before with order and orientation conservation

		if($n1!=$n2 && $$seen{$pair}!=2) { # the pair might still have conservation of order in another cog

		    while($j<$ncog1 && abs($n2-$n1)<=$genegap) {   # probably abs is not needed here due to the way I compute the id orders and also to the fact that the homology team program seems to return the witnesses in first chromosome in an increasing monotonical order of their start positions
			
			my $pair2=preserved_order($n1,$n2,$g1,$g2,$cog1[$i][1],$cog1[$j][1],\%cog3);
			if($pair2) {
			    if($$g1[$n1][1]==$$g1[$n2][1]) { # possible conserved Spair
				$ocsdpairs++;
				if(!$$seen{$pair}) { $csdpairs++;}
				if(abs($n2-$n1)==1) { 
				    $ocspairs++;
				    if(!$$seen{$pair}) { $cspairs++;}
				}
			    }
			    else {  # conserved pair of genes on different strands
				$ocddpairs++;
				if(!$$seen{$pair}) { $cddpairs++;}
				if(abs($n2-$n1)==1) { 
				    $ocdpairs++;
				    if(!$$seen{$pair}) { $cdpairs++;}
				}
			    }
			    $$seen{$pair}=2;
			    $$pairs{$pair}=$pair2;
			}
			elsif(!$$seen{$pair}) { # might still be conserved even if the order isn't
			    if($$g1[$n1][1]==$$g1[$n2][1]) { # possible conserved Spair
			    
				my $s2=$cog2{$cog1[$j][1]};
				if(!$s1 || !$s2 || $s2 == $s1) { # I have conserved direction
				    $csdpairs++;
				    if(abs($n2-$n1)==1) { $cspairs++;}
				    $$seen{$pair}=1;  # I MIGHT WANT TO ADD THE CONSERVED PAIR HERE
				}
			    
			    }
			    else {  # conserved pair of genes on different strands
				$cddpairs++;
				if(abs($n2-$n1)==1) { $cdpairs++;}
				$$seen{$pair}=1; # I MIGHT WANT TO ADD THE CONSERVED PAIR HERE
			    }
			}
			
			$j++;
			if($j<$ncog1) { 
			    $n2=$cog1[$j][0];
			    #$pair="$n1.$n2";
			    $pair=$$g1[$n1][3].".".$$g1[$n1][1].".".$$g1[$n2][3].".".$$g1[$n2][1];
			}
			
		    }
		}
	    }
	}
    }
    close(F);

    return($csdpairs,$cddpairs,$cspairs,$cdpairs,$ocsdpairs,$ocddpairs,$ocspairs,$ocdpairs);
}
			

sub compute_mult {
    my ($g,$circular,$bpgap,$genegap,$shift)=@_;

    my $ngenes=scalar(@{$g});

    my $ndir=1;
    my $sdpairs=0;
    my $ddpairs=0;
    my $spairs=0;
    my $dpairs=0;
    my $nstart=0;
    my $firstsign=0;
    my $lastsign;

    for(my $i=0;$i<$ngenes;$i++) {
	if($firstsign) { # seen first sign

	    if($$g[$i][0]-$$g[$i-1][2]>$bpgap) {
		$nstart=$i;
	    }
	    else {
		if($$g[$i][1]==$$g[$i-1][1]) {
		    $spairs++;
		}
		else {
		    $dpairs++;
		}
	    }
	    
	    if($$g[$i][1] != $lastsign) { 
		$ndir++;
	    }

	    my $max=$nstart > $i-$genegap ? $nstart : $i-$genegap;
	    for(my $j=$i-1;$j>=$max;$j--) {
		if($$g[$i][1]==$$g[$j][1]) {
		    $sdpairs++;
		}
		else {
		    $ddpairs++;
		}
	    }

	}
	else {
	    $firstsign=$$g[$i][1];
	}
	$lastsign=$$g[$i][1];

    }

    if($circular) {
	if($lastsign == $firstsign && $ndir>1) {
	    $ndir--;
	}

	if($$g[0][0]-$shift<=$bpgap) {
	    if($lastsign == $firstsign) { $spairs++;$sdpairs++;}
	    else { $dpairs++;$ddpairs++;}

	    my $j=$ngenes-1;
	    while($j>0 && $ngenes-$j<=$genegap && $$g[$j][0]-$$g[$j-1][2]<=$bpgap) {
		if($$g[0][1]==$$g[$j][1]) { $sdpairs++;}
		else { $ddpairs++;}
		$j--;
	    }

	    my $i=1;
	    while($i<$ngenes && $i<=$genegap-1 && $$g[$i][0]-$$g[$i-1][2]<=$bpgap) {
		my $gap=$i;
		$j=$ngenes-1;
		while($j>0 && $ngenes-$j+$gap<=$genegap && $$g[$j][0]-$$g[$j-1][2]<=$bpgap) {
		    if($$g[$i][1]==$$g[$j][1]) { $sdpairs++;}
		    else { $ddpairs++;}
		    $j--;
		}
		$i++;
	    }

	}
    }

    return($sdpairs,$ddpairs,$spairs,$dpairs,$ngenes,$ndir);

}

    


sub read_genome_cd {
    my ($pttfile,$g,$map,$id)=@_;

    my @k;
    my $lastgeneend=0;
    my $len;

    my $n=0;    
    my $m=0;
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

    my $end;

    while(<F>) {
	my @rec;
	my @a=split;
	my @b=split(/\./,$a[0]);

	my $keep=0;
	my $start=$b[0];
	$end=$b[-1];
	if($$map{$start}) { 
	    $start=$$map{$b[0]};
	    $keep=1;
	    $lastgeneend=$end;
	    $end=$$map{$b[-1]};
	}

	$rec[0]=$start;
	if($a[1] eq "+") { $rec[1]=1;}
	else { $rec[1]=-1};
	$rec[2]=$end;
	$rec[3]=$a[3];

	if($keep) { 
	    push(@{$k[$m]},@rec);
	    $m++;
	}
	else {
	    push(@{$$g[$n]},@rec);
	    $$id{$rec[0]}=$n;
	    $n++;
	}
    }
    close(F);

    for(my $j=0;$j<$m;$j++) {
	push(@{$$g[$n]},@{$k[$j]});
	$$id{$$g[$n][0]}=$n;
	$n++;
    }

    if($lastgeneend) { return($lastgeneend);}
    else { return($end-$len);}

}


sub readmap {
    my ($mapfile,$map)=@_;

    open(F,$mapfile);
    <F>;
    while(<F>) {
	chomp;
	my @a=split;
	$$map{$a[3]}=$a[1];
	$$map{$a[4]}=$a[2];
    }
    close(F);
}

sub preserved_order {
    my ($n11,$n12,$g1,$g2,$cogname1,$cogname2,$cog2)=@_;

    my $s11=$$g1[$n11][1];
    my $s12=$$g1[$n12][1];

    my $card1=scalar(@{$$cog2{$cogname1}});
    my $card2=scalar(@{$$cog2{$cogname2}});


    if($s11==$s12) {
	for(my $i=0;$i<$card1;$i++) {

	    my $n21=$$cog2{$cogname1}->[$i];
	    my $s21=$$g2[$n21][1];

	    for(my $j=0;$j<$card2;$j++) {

		my $n22=$$cog2{$cogname2}->[$j];
		my $s22=$$g2[$n22][1];
		
		if($s21==$s22 && ($n12-$n11)*($n22-$n21)*$s21*$s11>0) { 
		    my $pair2=$$g2[$n21][3]." ".$$g2[$n22][3];
		    return($pair2);
		}
	    }
	}
    }
    else {
	for(my $i=0;$i<$card1;$i++) {

	    my $n21=$$cog2{$cogname1}->[$i];
	    my $s21=$$g2[$n21][1];

	    for(my $j=0;$j<$card2;$j++) {

		my $n22=$$cog2{$cogname2}->[$j];
		my $s22=$$g2[$n22][1];
		
		if($s21!=$s22 && $s11*$s21*($n12-$n11)*($n22-$n21)>0) { 
		    my $pair2=$$g2[$n21][3]." ".$$g2[$n22][3];
		    return($pair2); 
		}
	    }

	}
    }
    return(0);
}
