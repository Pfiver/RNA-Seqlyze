#!/bin/bash -e

ADMIN_EMAIL=

TOPDIR=$PWD
PREFIX=/home/$USER/.local
WORKDIR=/home/$USER/rna-seqlyze-workdir
WORKDIR_DEV=$WORKDIR-dev

WWWDIR=/home/$USER/public_html
WWWBASE=/rna-seqlyze/
HOSTNAME=$(hostname -f)
GROUP=www-data
WORKER_USER=www-data
BIBODIR=$WORKDIR_DEV/buildbot
TRAC_DB=sqlite:///$WORKDIR_DEV/trac.db

# config variables documentation,
# copied from the 'Deployment' wiki page
doc_aADMIN_EMAIL='
   The email address of the application administrator;
   for example `admin@rna-seqlyze.com`.
'

doc_bTOPDIR='
   The root of the project source tree cloned from
   https://mu15ns00002.adm.ds.fhnw.ch/git/biocalc;
   for example `/home/user/rna-seqlyze`.
'
doc_cPREFIX='
   The project installation directory on the server;
   for example `/home/user/.local` or `/usr/local`, `/opt/biocalc`, etc.
'
doc_dWORKDIR='
   A directory on the server where lots of space should be available;
   for example`/home/user/rna-seqlyze-workdir`.
'
doc_eWORKDIR_DEV='
   A directory on the server where a second application instance
   along with a virtualenv in which it is run, is hosted;
   for example`/home/user/rna-seqlyze-workdir-dev`.
'

doc_fWWWDIR='
   The directory on the server containing the
   .htaccess file and the .wsgi scripts;
   for example `/home/user/public_html`.
'
doc_gWWWBASE='
   The path under which `WWWDIR` is accessible on the server
   ''from outside'' (e.g. with HTTP on port 80);
   for example `/rna-seqlyze/`.
'
doc_hHOSTNAME='
   The hostname under which the server is accessible from outside;
   for example `www.rna-seqlyze.com`.
'
doc_iGROUP='
   The unix group that the web application and the worker run as;
   for example `www-data`.
'
doc_jWORKER_USER='
   The unix user that the worker runs as;
   for example `www-data`.
'
doc_kTRAC_DB='
   A database url for a database where trac will keep its data;
   for example `sqlite:///$WORKDIR_DEV/trac.db`.
'
doc_lBIBODIR='
   A directory to hold one buildbot "master" and one buildbot
   "slave" base directory;
   for example `/home/user/buildbot`.
'

# os check
Debian= Ubuntu=
eval $(lsb_release -is)=true
if [ "$Debian$Ubuntu" != true ]
then
    echo Only Debian and Ubuntu OSs are supported so far.
    exit 1
fi

# dev check
if [ -n "$WORKDIR_DEV" ]
then
    echo "Do you want to install a development instance"
    read -p "in parallel with the production instance ? [Yes] " dev
    if ! [[ -z "$dev" -o "$dev" = [Yy]* ]]
    then
        unset WORKDIR_DEV TRAC_DB BIBODIR
        for doc in ${!doc_*}
        do
            var=${doc#doc_?}
            case $var in WORKDIR_DEV|TRAC_DB|BIBODIR)
                unset $doc
                break
                ;;
            esac
        done
    fi
fi

# if not manually configured at the top of the script,
# print the documentation for and ask for the value of each variable
if [ -z "$ADMIN_EMAIL" ]
then
    echo "You need to specify some configuration values."
    echo
    echo "To use the defaults, shown [in brackets], just hit return."
    for doc in ${!doc_*}
    do
        var=${doc#doc_?}
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

# WWWBASE: strip slashes
WWWBASE=${WWWBASE%/}
WWWBASE=${WWWBASE#/}

# helper functions
confsub=$(
    for doc in ${!doc_*}
    do
        var=${doc#doc_?}
        echo -n "s|@@$var@@|${!var//|/\\|}|g;"
    done
)
subcat() {
    sed "$confsub" "$1"
}
su() {
    if [ "$1" = -v ]
    then
        shift
        echo "$ ${Debian:+command su}${Ubuntu:+sudo sh} -c '$2'"
    fi
    ${Debian:+command su}${Ubuntu:+sudo sh} -c "$2"
}

# remember the current directory
CURDIR=$PWD
PATH=$PATH:$PREFIX/bin

# dependencies (su)
boost=libboost1.49-all-dev
[ -n "$(apt-cache search $boost)" ] || boost=
deps=(
    ntp ssh vim
    curl bzip2 sqlite3 deborphan python-lxml libapache2-mod-wsgi
    devscripts gcc-4.3-base
    cmake libbz2-dev libpng12-dev python-all-dev python3-all-dev $boost
)
su -v -c "apt-get install ${deps[*]}"

# boost 1.49 (su)
if [ -z "$boost" ]
then
    tmpdir=$(mktemp -d)
    base_url=http://ftp.ch.debian.org/debian/pool/main/b/boost1.49
    curl -JLOOO $base_url/boost1.49_1.49.0{.orig.tar.bz2,-3.1.{dsc,debian.tar.gz}}
    su -c "apt-get install libicu-dev mpi-default-dev bison flex \
                           docbook-to-man help2man xsltproc doxygen gccxml"
    cd $tmpdir
    dpkg-source -x boost1.49*.dsc
    cd boost1.49*
    dpkg-buildpackage -b
    cd ..
    su -v -c "dpkg -i *.deb"
    cd
    rm -rf $tmpdir
fi

# add $USER and $WORKER_USER to $GROUP (su)
gid=$(getent group $GROUP | cut -d : -f 3)
for user in $USER $WORKER_USER
do
    found=false
    for id in $(id -G)
    do
        if [ $id = $gid ]
        then
            found=true
            break
        fi
    done
    $found || su -v -c "usermod -a -G $GROUP $user"
done

# distribute
if ! python -c 'import setuptools' 2> /dev/null
then
    curl -O http://python-distribute.org/distribute_setup.py
    python << END_OF_PYTHON
import os, sys
os.makedirs(os.path.join('$PREFIX', 'lib',
    'python%d.%d' % sys.version_info[:2], 'site-packages'))
import distribute_setup as ds
ds._install(ds.download_setuptools(), ('--prefix', '$PREFIX'))
END_OF_PYTHON
    rm distribute_setup.py* distribute-*.tar.gz
fi

# docopt
if ! easy_install --prefix $PREFIX "docopt>0.4.1"
then
    tmpdir=$(mktemp -d)
    git clone https://github.com/docopt/docopt.git $tmpdir
    cd $tmpdir
    sed -i 's/version = "0.4.1/&.post1/' setup.py
    python setup.py install --prefix $PREFIX
    cd
    rm -rf $tmpdir
fi

# rna-seqlyze
cd $TOPDIR/src/rna-seqlyze
python setup.py install --prefix=$PREFIX

# 3rd-party software
cd $TOPDIR
rnas-install --prefix=$PREFIX

# apache/mod_wsgi environment
mkdir -p $WWWDIR
cd $TOPDIR/var/conf-files
subcat rna-seqlyze.wsgi > $WWWDIR/rna-seqlyze.wsgi
if [ -z "$WORKDIR_DEV" ]
then
    subcat htaccess > $WWWDIR/.htaccess
else
    subcat htaccess-dev >> $WWWDIR/.htaccess
    sed "s|@@WORKDIR@@|${WORKDIR_DEV//|/\\|}|;" \
        rna-seqlyze.wsgi > $WWWDIR/rna-seqlyze-dev.wsgi
fi

# utils
cd $TOPDIR/var/utils
for util in *
do
    subcat $util > $PREFIX/bin/$util
done

# crontab
cd $TOPDIR/var/conf-files
subcat crontab${WORKDIR_DEV:+-dev} | crontab

# workdirs
cd $TOPDIR/src/rna-seqlyze
rnas-init --group=$GROUP $WORKDIR
subcat rnaseqlyze.ini > $WORKDIR/rnaseqlyze.ini
if [ -n "$WORKDIR_DEV" ]
then
    rnas-init --group=$GROUP $WORKDIR_DEV
    ln -s $WORKDIR/shared_data $WORKDIR_DEV
    subcat rnaseqlyze.ini > $WORKDIR_DEV/rnaseqlyze.ini
fi

# trac
if [ -n "$TRAC_DB" ]
then
    # apache/mod_wsgi script
    cd $TOPDIR/var/conf-files
    for wsgi in {trac,debug}.wsgi
    do
        subcat $wsgi > $WWWDIR/$wsgi
    done
    ln -s trac.wsgi $WWWDIR/traclogin.wsgi

    # repo sync git hook
    cd $TOPDIR/.git
    for hook in post-commit post-merge
    do
        ln -vnfs $PREFIX/bin/trac-repo-sync hooks/$hook
    done

    # create trac.ini
    cd $TOPDIR/var/trac-env
    mkdir -p $WORKDIR-dev
    ti=$WORKDIR-dev/trac.ini
    subcat conf/trac.ini.tpl > $ti
    ln -s $ti conf/trac.ini
    chmod 664 $ti
    chgrp $GROUP $ti

    # set up user accounts
    hp=$WWWDIR/.htpasswd.trac
    touch $hp
    chmod 640 $hp
    chgrp $GROUP $hp
    echo Add some trac users.
    unset u
    while true
    do
        read -p "Username: " u
        [ -z "$u" ] && break
        while ! htpasswd $hp $u
        do
            echo try again
        done
        echo "added '$u'"
    done

    # initialize databse
    if [[ $TRAC_DB = sqlite:///* ]]
    then
        dbfile=${TRAC_DB#sqlite:///}
        sqlite3 $dbfile < $TOPDIR/var/trac-env/sqlite-db-backup.sql
        chmod 664 $dbfile
        chgrp $GROUP $dbfile
    fi

    # adjust log file permissions
    logfile=$WORKDIR-dev/trac.log
    touch $logfile
    chmod 664 $logfile
    chgrp $GROUP $logfile

    # adjust 'files' directory permissions
    mkdir -p files
    chmod 775 files
    chgrp $GROUP files
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
        buildslave create-slave buildslave localhost:9989 $HOSTNAME pass
    )
    cd $TOPDIR/var/conf-files
    subcat buildbot-master.cfg > $BIBODIR/buildmaster/master.cfg
fi

# apache and -worker daemon (su)
cd $CURDIR
a2conf=rna-seqlyze-a2conf${WORKDIR_DEV:+-dev}
subcat $TOPDIR/var/conf-files/$a2conf > rna-seqlyze-a2conf
apache_conf='
chown root. rna-seqlyze-a2conf
mv rna-seqlyze-a2conf /etc/apache2/conf.d/rna-seqlyze
a2enmod proxy rewrite wsgi
service apache2 restart
'
subcat $TOPDIR/src/rna-seqlyze-worker/rna-seqlyze-service > rna-seqlyze-service
service_conf='
chown root. rna-seqlyze-service
chmod 755 rna-seqlyze-service
mv rna-seqlyze-service /etc/init.d/rna-seqlyze
insserv rna-seqlyze
service rna-seqlyze start
'
cat << END_OF_ROOT

Almost done!

There are two more things left to be done
that need to be carried out with root permissions:

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
if ! su -c "$apache_conf$service_conf"
then
    echo "$apache_conf$service_conf" > finish-installation
    echo "An error happend trying to issue the mentioned commands."
    echo "They have been written to a file called 'finish-installation'."
fi

# done
cat << END_OF_DONE

 RNA-Seqlyze installed in $PREFIX

 WORKDIR=$WORKDIR${WORKDIR_DEV:+
 WORKDIR_DEV=$WORKDIR_DEV}

END_OF_DONE
