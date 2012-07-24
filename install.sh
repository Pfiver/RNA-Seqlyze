#!/bin/bash -e

TOPDIR=/home/user/rna-seqlyze
PREFIX=/home/user/.local
WORKDIR=/home/user/rna-seqlyze-workdir
WWWDIR=/home/user/public_html
WWWBASE=/rna-seqlyze/

# optional
#BIBODIR=/home/user/buildbot

if [ $TOPDIR = /home/user/rna-seqlyze ]; then cat << 'END_OF_DOC'

To install rna-seqlyze, edit the variables on top of this script:

 `TOPDIR`:
   The root of the project source tree cloned from
   https://mu15ns00002.adm.ds.fhnw.ch/git/biocalc;
   for example `/home/user/rna-seqlyze`.

 `PREFIX`:
   The project installation directory on the server;
   for example `/home/user/.local` or `/usr/local`, `/opt/biocalc`, etc.

 `WORKDIR`:
   A directory on the server where lots of space should be available;
   for example`/home/user/rna-seqlyze-workdir`.

 `WWWDIR`:
   The directory on the server containing the
   .htaccess file and the .wsgi scripts;
   for example `/home/user/public_html`.

 `WWWBASE`:
   The path under which `WWWDIR` is accessible on the server
   ''from outside'' (e.g. with HTTP on port 80);
   for example `/rna-seqlyze/`.

 `BIBODIR`:
   A directory to hold one buildbot "master" and one buildbot
   "slave" base directory;
   for example `/home/user/buildbot`.

(copied from trac/wiki/Deployment)

END_OF_DOC
exit 1
fi

# helper function
install()
{
    sed "s|@@TOPDIR@@|$TOPDIR|;
         s|@@PREFIX@@|$PREFIX|;
         s|@@WORKDIR@@|$WORKDIR|;
         s|@@WWWDIR@@|$WWWDIR|;
         s|@@WWWBASE@@|$WWWBASE|;
         s|@@BIBODIR@@|$BIBODIR|;" $1 > $2
}

# apache/mod_wsgi environment
mkdir -p $WWWDIR
cd $TOPDIR/var/conf-files
install htaccess $WWWDIR/.htaccess
cd $TOPDIR/var/utils
for wsgi in *.wsgi
do
    install $wsgi $WWWDIR/$wsgi
done
sed "s|@@WORKDIR@@|$WORKDIR-dev|;" \
    rna-seqlyze.wsgi > $WWWDIR/rna-seqlyze-dev.wsgi

# utils
cd $TOPDIR/var/utils
shopt -s extglob
for bin in !(*.wsgi)
do
    install $bin $PREFIX/bin/$bin
done

if [ -n "$BIBODIR" ]
then
    # buildbot config
    cd $TOPDIR/var/conf-files
    install buildbot-master.cfg $BIBODIR/buildmaster/master.cfg
fi

# crontab
cd $TOPDIR/var/conf-files
crontab crontab 

# git hooks
cd $TOPDIR/.git
for hook in post-commit post-merge
do
    ln -vnfs $PREFIX/bin/trac-repo-sync hooks/$hook
done
