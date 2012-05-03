#!/bin/sh

if false
then
	echo Content-Type: text/plain
	echo
	echo $0
	echo $@
	env
	id
	echo
	RU=${REQUEST_URI%%\?*}
	RU=${RU//%??/X}
	echo $RU
	echo ${RU::${#RU}-${#PATH_INFO}+1}
	exit 0
fi

RU=${REQUEST_URI%%\?*}
RU=${RU//%??/X}
SCRIPT_NAME=${RU::${#RU}-${#PATH_INFO}+1}

export REMOTE_USER=$REDIRECT_REMOTE_USER

export PYTHONPATH=/home/aeshoh6c/lib/python2.7/site-packages

exec /home/aeshoh6c/bin/python2.7 cgi-bin/trac.cgi
