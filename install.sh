#!/bin/bash -e

TOPDIR=/home/$USER/rna-seqlyze
PREFIX=/home/$USER/.local
WORKDIR=/home/$USER/rna-seqlyze-workdir
WWWDIR=/home/$USER/public_html
WWWBASE=/rna-seqlyze/
BIBODIR=/home/$USER/buildbot

HOSTNAME=$(hostname -f)
GROUP=www-data
WORKER_USER=www-data
ADMIN_EMAIL=
TRAC_DB=sqlite:///$WORKDIR-dev/trac.db

# config variables documentation,
# copied from the 'Deployment' wiki page
doc_TOPDIR='
   The root of the project source tree cloned from
   https://mu15ns00002.adm.ds.fhnw.ch/git/biocalc;
   for example `/home/user/rna-seqlyze`.
'
doc_PREFIX='
   The project installation directory on the server;
   for example `/home/user/.local` or `/usr/local`, `/opt/biocalc`, etc.
'
doc_WORKDIR='
   A directory on the server where lots of space should be available;
   for example`/home/user/rna-seqlyze-workdir`.
'
doc_WWWDIR='
   The directory on the server containing the
   .htaccess file and the .wsgi scripts;
   for example `/home/user/public_html`.
'
doc_WWWBASE='
   The path under which `WWWDIR` is accessible on the server
   ''from outside'' (e.g. with HTTP on port 80);
   for example `/rna-seqlyze/`.
'
doc_BIBODIR='
   A directory to hold one buildbot "master" and one buildbot
   "slave" base directory;
   for example `/home/user/buildbot`.
'

doc_HOSTNAME='
   The hostname under which the server is accessible from outside;
   for example `www.rna-seqlyze.com`.
'
doc_GROUP='
   The unix group that the web application and the worker run as;
   for example `www-data`.
'
doc_WORKER_USER='
   The unix user that the worker runs as;
   for example `www-data`.
'
doc_ADMIN_EMAIL='
   The email address of the application administrator;
   for example `admin@rna-seqlyze.com`.
'
doc_TRAC_DB='
   A database url for a database where trac will keep its data;
   for example `sqlite:///$WORKDIR-dev/trac.db`.
'

# os check
Debian= Ubuntu=
eval $(lsb_release -is)=true
if [ "$Debian$Ubuntu" != true ]
then
    echo Only Debian and Ubuntu OSs are supported so far.
    exit 1
fi

# if not already manually configured at the top of the script,
# print the documentation for and ask for the value of each variable
if [ -z "$ADMIN_EMAIL" ]
then
    echo "You need to specify some configuration values."
    echo
    echo "To use the defaults, shown [in brackets], just hit return."
    for doc in ${!doc_*}
    do
        var=${doc#doc_}
        echo
        echo "${!doc}"
        while true
        do
            read -p "$var = [${!var}] " val || exit
            [ -n "$val" ] && eval $var=\$val
            [ -n "${!var}" ] && break
            echo "$var can't be empty."
        done
    done
    echo
fi

# helper function
confsub=$(
    for doc in ${!doc_*}
    do
        var=${doc#doc_}
        echo -n "s|@@$var@@|${!var//|/\\|}|g;"
    done
)
subcat() { sed "$confsub" "$1"; }

# remember the current directory
CURDIR=$PWD

# rna-seqlyze
cd $TOPDIR/src/rna-seqlyze
python setup.py install --prefix=$PREFIX

# 3rd-party software
cd $TOPDIR
rnas-install --prefix=$PREFIX

# apache/mod_wsgi environment
mkdir -p $WWWDIR
cd $TOPDIR/var/conf-files
subcat htaccess > $WWWDIR/.htaccess
for wsgi in *.wsgi
do
    subcat $wsgi > $WWWDIR/$wsgi
done
if [ -n "$TRAC_DB" ]
then
    ln -s trac.wsgi $WWWDIR/traclogin.wsgi
    sed "s|@@WORKDIR@@|$WORKDIR-dev|;" \
        rna-seqlyze.wsgi > $WWWDIR/rna-seqlyze-dev.wsgi
fi

# utils
cd $TOPDIR/var/utils
for util in *
do
    subcat $util > $PREFIX/bin/$util
done

# buildbot
if [ -n "$BIBODIR" ]
then
    cd $TOPDIR/var/conf-files
    subcat buildbot-master.cfg > $BIBODIR/buildmaster/master.cfg
fi

# crontab
{
    cat << END_OF_CRONTAB 
SHELL=/bin/bash

rnas_workdir=$WORKDIR
rnas_bibodir=$BIBODIR
rnas_wwwbase=$WWWBASE

# m  h   dom mon dow command
 12  *   *   *   *   . $TOPDIR/bash-env; update-available-list
END_OF_CRONTAB
    if [ -n "$TRAC_DB" ]
    then
        echo " 13  *   *   *   *  " \
                "rnas_workdir=$WORKDIR-dev . $TOPDIR/bash-env; cron-hourly"
    fi
} | crontab 

# git hooks
cd $TOPDIR/.git
for hook in post-commit post-merge
do
    ln -vnfs $PREFIX/bin/trac-repo-sync hooks/$hook
done

# workdirs
rnas-init $WORKDIR
rnas-init $WORKDIR-dev
cd $TOPDIR/src/rna-seqlyze
subcat rnaseqlyze.ini > $WORKDIR/rnaseqlyze.ini
subcat rnaseqlyze.ini > $WORKDIR-dev/rnaseqlyze.ini

# trac
if [ -n "$TRAC_DB" ]
then
    # create trac.ini
    cd $PREFIX/var/trac-env
    mkdir -p $WORKDIR-dev
    ti=$WORKDIR-dev/trac.ini
    subcat conf/trac.ini.tpl > $ti
    ln -s $ti conf/trac.ini
    chmod 664 $ti
    chgrp $GROUP $ti

    # adjust 'files' directory permissions
    mkdir -p files
    chmod 775 files
    chgrp $GROUP files

    # set up user accounts
    hp=$WWWDIR/.htpasswd.trac
    touch $hp
    echo Add some trac users.
    unset u
    while true
    do
        read -p "Username: " u
        [ -z "$u" ] && break
        htpasswd -b $hp $u
    done

    # initialize databse
    if [[ $TRAC_DB = sqlite:///* ]]
    then
        sqlite3 ${TRAC_DB#sqlite:///} < \
            $TOPDIR/var/trac-env/sqlite-db-backup.sql
    fi
fi

# buildbot
#  http://pypi.python.org/pypi/virtualenv
#  http://buildbot.net/buildbot/docs/current/manual/installation.html
if [ -n "$BIBODIR" ]
then
    easy_install --user virtualenv
    virtualenv --system-site-packages --distribute $BIBODIR
    (
        cd $BIBODIR
        . bin/activate
        pip install buildbot
        pip install buildbot-slave
        buildbot create-master -r buildmaster
        buildslave create-slave buildslave localhost:9989 biopython-slave pass
    )
fi

# apache, worker service
cd $CURDIR
subcat $TOPDIR/src/rna-seqlyze-worker/rna-seqlyze-service > rna-seqlyze-service
subcat $TOPDIR/var/conf-files/rna-seqlyze-a2conf > rna-seqlyze-a2conf

# as root
apache_conf='
mv rna-seqlyze-a2conf /etc/apache2/conf.d/rna-seqlyze
a2enmod proxy rewrite wsgi
service apache2 restart
'
service_conf='
mv rna-seqlyze-service /etc/init.d/rna-seqlyze
insserv rna-seqlyze
service rna-seqlyze start
'
cat << END_OF_ROOT

Almost done!

There are two things left to be done
that need to be carried out with _root permissions_:

1. To configure the apache webserver, the file 'rna-seqlyze-a2conf',
   that has just been created in the current directory, must be moved to
   '/etc/apache2/conf.d'. Then all required modules must be activated and
   the apache daemon needs to be restarted by issuing the following commands:

----
$apache_conf
----

2. To configure the rna-seqlyze-worker daemon (service), the file
   'rna-seqlyze-service', that has just been created in the current
   directory, must be moved to '/etc/init.d'. Then the service must be
   activated and started by issuing the following commands:

----
$service_conf
----

Enter ${Debian:+the root}${Ubuntu:+your} password to issue those commands now.

END_OF_ROOT
if ! eval ${Debian:+su -c}${Ubuntu:+sudo sh -c} "$apache_conf$service_conf"
then
    echo "$apache_conf$service_conf" > finish-installation
    echo "An error happend trying to issue the mentioned commands."
    echo "They have been written to a file called 'finish-installation'."
fi

# done
cat << END_OF_DONE

 RNA-Seqlyze installed in $PREFIX

 WORKDIR=$WORKDIR

END_OF_DONE
