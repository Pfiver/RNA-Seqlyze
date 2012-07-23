#!/bin/bash

echo installing RNA-Seqlyze....
echo

echo ...implementation comming soon
exit

TOPDIR=
PREFIX=
WORKDIR=
WORKDIR_DEV=

WEBDIR=

confsub=""
confsub+=" s|@@TOPDIR@@|$TOPDIR|;"
confsub+=" s|@@PREFIX@@|$PREFIX|;"
confsub+=" s|@@WORKDIR@@|$WORKDIR|;"
confsub+=" s|@@WORKDIR_DEV@@|$WORKDIR_DEV|;"
confsub+=" s|@@WEBDIR@@|$WEBDIR|;"

cd ~

if $git_init
	git init
	echo $TOPDIR/.git/objects > .git/objects/info/alternates
	git read-tree $(cd $TOPDIR; git rev-parse HEAD:var/biocalc-env)
	git checkout .
	sed -i $confsub $(git ls-files)
else
	echo "(not overwriting existing files)"
	git --git-dir $TOPDIR/.git archive HEAD:var/biocalc-env | tar kvx
	sed -i $confsub $(git --git-dir $TOPDIR/.git \
				ls-tree --name-only -r HEAD:var/biocalc-env)
fi
