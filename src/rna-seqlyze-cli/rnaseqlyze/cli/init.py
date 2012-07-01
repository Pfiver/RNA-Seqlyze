"""\
RNA-Seqlyze Init

(Re-)initialize an rnaseqlyze 'workdir'.

Usage:
    rnas-init [--recreatedb] <workdir>
    rnas-init --development <workdir>
    rnas-init -h|--help

Options:
    --recreatedb    Remove and re-initialize the database if it exists.

    --development   Use the development versions of the config file templates.

Arguments:
    <workdir>       The filesystem path to the directory to be initialized.
                    If the directory already exists, by default, existing
                    files inside that directory are not overwritten.

Documentation:

    The 'workdir' holds

        - configuration files           *.ini
        - the application database      rnaseqlyze.db
        - log files                     *.log
        - shared data                   shared_data/
        - individual analysis data      analyses/


    The rnas-init command

        - creates the <workdir> if it does not already exist

        - copies default configuration files 'rnaseqlyze.ini', 'web.ini'
          and 'worker.ini' into the <workdir> if they do not already exist.

        - initializes the database that is configured in the 'rnaseqlyze.ini'
          config file if it doesn't already exist and --recreatedb is not given.


    The database to be initialized is configured with "db_url =" in the
    "[rnaseqlyze]" section in 'rnaseqlyze.ini'.  It is expected to be an sqlite
    database.

    If the command creates the sqlite database file, it changes it's unix
    access mode to (octal) 0664 and the group membership is changed to <group>.
    <group> can be confgured in 'rnaseqlyze.ini'.

    If the command creates the <workdir>, it changes it's unix access mode to
    (octal) 0775 and the group membership is also changed to <group>.
"""

from __future__ import print_function

import os, sys, grp, shutil

import pkg_resources
from sqlalchemy import create_engine

import rnaseqlyze
import rnaseqlyze.web
import rnaseqlyze.worker
from rnaseqlyze.core.orm import Entity

def main():

    import docopt
    opts = docopt.docopt(__doc__)

    workdir = opts['<workdir>']

    # create the workdir if it does not exist
    wd_created = False
    if not os.path.isdir(workdir):
        if os.path.exists(workdir):
            print("not a directory: '%s'" % workdir, file=sys.stderr)
            sys.exit(1)
        os.makedirs(workdir)
        wd_created = True

    # create each config file that does not exist
    for pkg in rnaseqlyze, rnaseqlyze.web, rnaseqlyze.worker:

        conf_name = pkg.__name__.split('.')[-1] + ".ini"
        conf_path = os.path.join(wordir, conf_name)
        if os.path.exists(conf_path):
            continue

        ini = opts['--development'] and "development.ini" or "production.ini"
        res = pkg_resources.resource_stream(pkg.project_name, ini)
        shutil.copyfileobj(res, open(conf_path, "w"))

    # init rnaseqlyze configuration
    rnaseqlyze.configure(os.path.join(workdir, 'rnaseqlyze.ini'))

    # delayed because 'rnaseqlyze.group' was
    # not known before calling rnaseqlyze.configure() above
    if wd_created:
        # change permission bits
        os.chmod(workdir, 0775)
        # change group membership
        os.chown(workdir, -1, grp.getgrnam(rnaseqlyze.group).gr_gid)

    # get the database file path
    db_path = rnaseqlyze.db_url.split(":", 1)[1]

    # remove the databse file
    # if it exists and --recreatedb is given
    if os.path.exists(db_path) and opts['--recreatedb']:
        os.unlink(db_path)

    # create the database if it doesn't exist
    if not os.path.exists(db_path):

        # create the file and initialize the schema
        Entity.metadata.create_all(create_engine(rnaseqlyze.db_url))

        # change permission bits
        os.chmod(db_path, 0664)
        # change group membership
        os.chown(db_path, -1, grp.getgrnam(rnaseqlyze.group).gr_gid)

    print "workdir initialized"
