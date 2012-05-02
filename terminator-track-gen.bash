#!/bin/bash

gb=$1

id=${gb##*/} id=${id%.gb}

tmpdir=/tmp/terminator-track-gen

[ -d "$tmpdir" ] || mkdir $tmpdir

t=$tmpdir/$id.terminators
fa=$tmpdir/$id.fasta
ptt=$tmpdir/$id.ptt

set -C # If set, bash does not overwrite an existing
       # file with the >, >&, and <> redirection operators

gb2ptt.py $gb > $ptt
gb2fasta.py $gb > $fa
sed -i "1s/.*/>$id/" $fa
transterm.bash --all-context $fa $ptt > $t

        # goal:
        #       track name="BED" description="BED Test Track" visibility=2 itemRgb="on" 
        #       chr     1       100     Foo     500     +       10      90      255,0,0
        #       chr     200     250     Bar     900     +       205     245     0,255,0
        #       chr     400     425     Baz     800     -       405     420     0,0,255

{
	echo 'track name="BED" description="BED Test Track" visibility=2 itemRgb="on"'

	< $t grep "^  TERM" |
		while read line
	do
		id=${line:5:6} str=${line:31:1} con=${line:37:3} col=$((con*128))
		if [ "$str" = "+" ]
		then beg=${line:12:7} end=${line:22:7}
		else beg=${line:22:7} end=${line:12:7}
		fi
		echo "chr	$beg	$end	TERM_$id	$con	$str	$beg	$end	$col,$col,$col"
	done
} \
	> $id-transterm-all.bed
