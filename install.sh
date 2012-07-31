#!/bin/bash -e

# variable defaults
#  if ADMIN_EMAIL is set as well the script will not ask for any values
#  if WORKDIR_DEV is not set, no development instance will be installed
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
WORKER_PORT=5433
WORKER_PORT_DEV=6544
BIBODIR=$WORKDIR_DEV/buildbot
TRAC_DB=sqlite:///$WORKDIR_DEV/trac.db

# variables documentation
#  copied from 'Deployment' wiki page
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
doc_kWORKER_PORT='
   The tcp listening port of the worker daemon
'
doc_lWORKER_PORT_DEV='
   The tcp listening port of the development worker daemon
'
doc_mTRAC_DB='
   A database url for a database where trac will keep its data;
   for example `sqlite:///$WORKDIR_DEV/trac.db`.
'
doc_nBIBODIR='
   A directory to hold one buildbot "master" and one buildbot
   "slave" base directory;
   for example `/home/user/buildbot`.
'

# remember the current directory
CURDIR=$PWD

# os check
Debian= Ubuntu=
eval $(lsb_release -is)=true
if [ "$Debian$Ubuntu" != true ]
then
    echo Only Debian and Ubuntu OSs are supported so far.
    exit 1
fi

# 'python check'
PYVER=$(python -c 'import sys; print sys.version[:3]')

# if not manually configured at the top of the script,
# print the documentation for and ask for the value of each variable
if [ -z "$ADMIN_EMAIL" ]
then
    # dev check
    if [ -n "$WORKDIR_DEV" ]
    then
        cat << 'END_OF_QUESTION'

  Do you want to install a development instance
  in parallel with the production instance ?
  If you answer yes, a 'trac' and a 'buildbot' instance will
  be configured and installed under $WORKDIR_DEV as well.

END_OF_QUESTION
        read -p "Yes/No [Yes] ? " dev
        if ! [[ -z "$dev" || "$dev" = [Yy]* ]]
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

    # ask for config values
    cat << 'END_OF_MESSAGE'

  You need to specify some configuration values.
  To use the defaults, shown [in brackets], just hit return.
END_OF_MESSAGE
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

# install dev env ?
devinst=${WORKDIR_DEV:+true}
: ${devinst:=false}

# WWWBASE: strip slashes
WWWBASE=${WWWBASE%/}
WWWBASE=${WWWBASE#/}

# confsub script
confsub=$(mktemp)
chmod 755 $confsub
echo '#!/bin/sed -f' > $confsub
for doc in ${!doc_*}
do
    var=${doc#doc_?}
    echo "s|@@$var@@|${!var//|/\\|}|g"
done \
    >> $confsub

# use binaries in $PREFIX/bin
PATH=$PREFIX/bin:$PATH

# use modules in and create directory
# $PREFIX/lib/python$PYVER/lib/site-packages
PYSITE=lib/python$PYVER/site-packages
PYTHONPATH=$PREFIX/$PYSITE
mkdir -p $PYTHONPATH

# @@PYTHON_PATH@@ confsub
#  required in rna-seqlyze-a2conf*
echo "s|@@PYTHON_PATH@@|${PYTHONPATH//|/\\|}|" >> $confsub

# devinst PYTHON_PATH
if $devinst
then
    PYTHONPATH_DEV=$WORKDIR_DEV/$PYSITE:$PYTHONPATH
    echo "s|@@PYTHON_PATH_DEV@@|${PYTHONPATH_DEV//|/\\|}|" >> $confsub
fi

# su/sudo helper
su() {
    if [ "$1" = -v ]
    then
        shift
        echo "$ ${Debian:+command su}${Ubuntu:+sudo sh} -c '$2'"
    fi
    ${Debian:+command su}${Ubuntu:+sudo sh} -c "$2"
}

# dependencies (su)
boost=libboost1.49-all-dev
[ -n "$(apt-cache search $boost)" ] || boost=
deps=(
    ntp ssh vim
    curl bzip2 sqlite3 deborphan python-lxml libapache2-mod-wsgi
    devscripts gcc-4.3-base gccxml
    cmake libbz2-dev libpng12-dev python-all-dev python3-all-dev $boost
)
echo root permissions required to install dependenceies:
su -v -c "apt-get install ${deps[*]}"

# boost 1.49 (su)
if [ -z "$boost" ]
then
    tmpdir=$(mktemp -d)
    base=http://ftp.ch.debian.org/debian/pool/main/b/boost1.49
    curl -JLOOO $base/boost1.49_1.49.0{.orig.tar.bz2,-3.1.{dsc,debian.tar.gz}}
    echo "root permissions required to install 'boost' 1.49 build dependencies"
    su -v -c "apt-get install libicu-dev mpi-default-dev bison flex \
                              docbook-to-man help2man xsltproc doxygen"
    cd $tmpdir
    dpkg-source -x boost1.49*.dsc
    cd boost1.49*
    dpkg-buildpackage -b
    cd ..
    echo "root permissions required to install locally built 'boost' 1.49"
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
    if ! $found
    then
        echo "root permissions required to add '$user' to '$GROUP' group"
        su -v -c "usermod -a -G $GROUP $user"
    fi
done

# distribute
if ! python -c 'import setuptools' 2> /dev/null
then
    curl -O http://python-distribute.org/distribute_setup.py
    python << END_OF_PYTHON
import distribute_setup as ds
ds._install(ds.download_setuptools(), ('--prefix', '$PREFIX'))
END_OF_PYTHON
    rm distribute_setup.py* distribute-*.tar.gz
fi

# docopt
if ! easy_install --prefix $PREFIX "docopt > 0.4.1"
then
    tmpdir=$(mktemp -d)
    git clone https://github.com/docopt/docopt.git $tmpdir
    cd $tmpdir
    sed -i 's/version = "0.4.1/&.post1/' setup.py
    python setup.py install --prefix $PREFIX
    cd
    rm -rf $tmpdir
fi

if $devinst
then
    # virtualenv
    easy_install --prefix $PREFIX virtualenv

    # GitPython
    #  - git-ls-authors dependency
    #  - need a version that includes commit 864cf1a4
    #    (https://github.com/gitpython-developers/GitPython/commit/864cf1a4)
    if ! python << 'END_OF_PYTHON'
import sys, pkg_resources
try: pkg_resources.require('GitPython > 0.3.1')
except: sys.exit(1)
END_OF_PYTHON
    then
        tmpdir=$(mktemp -d)
        git clone https://github.com/gitpython-developers/GitPython.git $tmpdir
        cd $tmpdir
        git checkout 864cf1a4
        echo 0.3.1.post1 > VERSION
        rm -rf git/test/fixtures
        python setup.py install --prefix $PREFIX
        cd
        rm -rf $tmpdir
    fi
fi

# rna-seqlyze
for dir in $TOPDIR/src/rna-seqlyze*
do
    (cd $dir && python setup.py install --prefix $PREFIX)
done

# 3rd-party software
cd $TOPDIR
rnas-setup --prefix $PREFIX

# apache/mod_wsgi environment
mkdir -p $WWWDIR
cd $TOPDIR/var/conf-files
$confsub rna-seqlyze.wsgi > $WWWDIR/rna-seqlyze.wsgi
if ! $devinst
then
    $confsub htaccess > $WWWDIR/.htaccess
else
    $confsub htaccess-dev >> $WWWDIR/.htaccess
    sed "s|@@WORKDIR@@|@@WORKDIR_DEV@@|;" rna-seqlyze.wsgi |
        $confsub > $WWWDIR/rna-seqlyze-dev.wsgi
fi

# utils
cd $TOPDIR/var/utils
for util in *
do
    $confsub $util > $PREFIX/bin/$util
    chmod 755 $PREFIX/bin/$util
done

# crontab
if $devinst
then
    cd $TOPDIR/var/conf-files
    {
        crontab -l
        echo
        echo \# --- added by $0 on $(date) ---
        $confsub crontab
    } |
        crontab
fi

# workdirs
cd $TOPDIR/src/rna-seqlyze
rnas-init --group=$GROUP $WORKDIR
$confsub rnaseqlyze.ini > $WORKDIR/rnaseqlyze.ini
if $devinst
then
    rnas-init --group=$GROUP \
        --development $WORKDIR_DEV
    ln -s $WORKDIR/shared_data $WORKDIR_DEV
    sed 's|@@WORKER_PORT@@|@@WORKER_PORT_DEV@@|g' rnaseqlyze.ini |
        $confsub > $WORKDIR_DEV/rnaseqlyze.ini

    cd $TOPDIR/var/conf-files
    $confsub bash-env > $WORKDIR_DEV/bash-env

    # devinst virtualenv
    virtualenv --system-site-packages --distribute $WORKDIR_DEV
# doesn't work as intended
#    (
#        . $WORKDIR_DEV/bin/activate
#        for dir in $TOPDIR/src/rna-seqlyze*
#        do
#            (cd $dir && python setup.py develop)
#        done
#    )
# doing it by hand instead
    cd $WORKDIR_DEV/lib/python$PYVER/site-packages
    mv easy-install.pth{,~}
    {
        head -n-1 easy-install.pth~
        ls -d1 $TOPDIR/src/rna-seqlyze*
        tail -n-1 easy-install.pth~
    } \
        > easy-install.pth
    rm easy-install.pth~
fi

# trac
if $devinst
then
    # apache/mod_wsgi script
    cd $TOPDIR/var/conf-files
    for wsgi in {trac,debug}.wsgi
    do
        $confsub $wsgi > $WWWDIR/$wsgi
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
    mkdir -p $WORKDIR_DEV
    ti=$WORKDIR_DEV/trac.ini
    $confsub conf/trac.ini.tpl > $ti
    ln -vnsf $ti conf/trac.ini
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
    logfile=$WORKDIR_DEV/trac.log
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
if $devinst
then
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
    $confsub buildbot-master.cfg > $BIBODIR/buildmaster/master.cfg
fi

# apache and -worker daemon (su)
cd $CURDIR
a2conf=rna-seqlyze-a2conf${WORKDIR_DEV:+-dev}
$confsub $TOPDIR/var/conf-files/$a2conf > rna-seqlyze-a2conf
apache_conf='
mv rna-seqlyze-a2conf /etc/apache2/conf.d/rna-seqlyze
chown root. /etc/apache2/conf.d/rna-seqlyze
a2enmod proxy proxy_http rewrite wsgi
service apache2 restart
'
$confsub $TOPDIR/var/conf-files/rna-seqlyze-service > rna-seqlyze-service
service_conf='
mv rna-seqlyze-service /etc/init.d/rna-seqlyze
chown root. /etc/init.d/rna-seqlyze
chmod 755 /etc/init.d/rna-seqlyze
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

# save confsub script
if $devinst
then
    mv $confsub $WORKDIR_DEV/confsub
else
    rm $confsub
fi

# done
cat << END_OF_DONE

 RNA-Seqlyze installed in $PREFIX

 WORKDIR=$WORKDIR${WORKDIR_DEV:+
 WORKDIR_DEV=$WORKDIR_DEV}

END_OF_DONE
