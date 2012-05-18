#!/bin/bash

all=
if [[ "$1" = -a* ]]
then
	shift
	all=\ --all-context
fi

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
transterm.bash$all $fa $ptt > $t

        # goal:
        #       track name="BED" description="BED Test Track" visibility=2 itemRgb="on" 
        #       chr     1       100     Foo     500     +       10      90      255,0,0
        #       chr     200     250     Bar     900     +       205     245     0,255,0
        #       chr     400     425     Baz     800     -       405     420     0,0,255

{
	echo 'track name="TransTermHP'"$all"'" description="TransTermHP Terminator Predictions'"$all"'" visibility=2 itemRgb="on"'

	< $t grep "^  TERM" |
		while read term id beg dash end str pos con rest
	do
		col=$((100-con))	# let the color vary between 0 (black) and 100 (gray)
		if [ "$str" = "-" ]	# switch begin & end if on reverse strand
		then tmp=$beg beg=$end end=$tmp
		fi
		echo "chr	$beg	$end	TERM_$id	$con	$str	$beg	$end	$col,$col,$col"
	done
}

echo "temporary files stored in '$tmpdir'" >&2
