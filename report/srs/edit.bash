#!/bin/bash

if ! [ -x $(type -p inotifywait) ]
then
	echo install inotify-tools: $ sudo apt-get install inotify-tools
	echo "> gives you a directory watcher command line ulity"
	exit
fi

if ! make					# compile once
then
	echo install some latex stuff: $ sudo apt-get install \
		texlive-latex-recommended texlive-humanities texlive-latex-extra
	echo "> gives you koma-script, lineno, lastpage & todonotes"
	exit
fi

						# determine version
VERSION=$(make -p . | sed '/^VERSION/!d;s/.*= //')
pdf=SRS-$VERSION.pdf
evince "$pdf" &					# start evince & gedit

inputs=(*.tex *.cls)
gedit "${inputs[@]}" &

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
