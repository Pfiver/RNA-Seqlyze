#!/bin/bash

if ! [ -x $(type -p inotifywait) ]
then
        echo
	echo "  install inotify-tools:"
        echo
        echo " $ sudo apt-get install inotify-tools"
        echo
	echo "      -> gives you a directory watcher command line ulity"
	exit
fi

if ! make					# compile once
then
        echo
	echo "  install some latex stuff:"
        echo
        echo " $ sudo apt-get install texlive-latex-recommended \\"
        echo "                         texlive-humanities texlive-latex-extra"
        echo
	echo "      -> gives you koma-script, lineno, lastpage & todonotes"
	exit
fi

start=true
case $1 in ?*) start=false;; esac
                                                # determine version
VERSION=$(make -p . | sed '/^VERSION/!d;s/.*= //')
pdf=BachelorThesis-$VERSION.pdf
$start && evince "$pdf" &                       # start evince & gedit

inputs=(*.tex *.cls)
$start && gedit "${inputs[@]}" &

exec 10< <(inotifywait -m --format "%e %f" .)	# watch the current directory
while read -u 10 event arg			# read events from inotifywait
do
#	echo $event $arg
	for i in "${inputs[@]}"			# check if the event affects US
	do
		if [ "$i" = "$arg" ]
		then
			IFS=,
			for e in $event		# check if it's a "modify" event
			do
				if [ "$e" = MODIFY ] \
				   || [ "$e" = MOVED_TO ]
				then
					make	# recompile the pdf from sources
				fi
			done
			IFS=$' \t\n'
		fi
	done
done
