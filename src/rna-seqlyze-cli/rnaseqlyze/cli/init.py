"""\
RNA-Seqlyze Init

(Re-)initialize an rnaseqlyze 'workdir'.

Usage:
    rnas-init <workdir>
    rnas-init --recreatedb <workdir>
    rnas-init --development <workdir>
    rnas-init -h|--help

Options:
    --recreatedb    Remove and re-initialize the database if it exists.

    --development   Use the development versions of the config file templates.

Arguments:
    <workdir>       The filesystem path to the directory to be initialized.
                    If the directory already exists, by default, existing
                    files inside that directory are not overwritten.

                 + ---------------------------------------------------------- +
                 |  The `WORKDIR` variable in the "/etc/init.d/rnaseqlyze.sh" |
                 |  worker daemon startup script and the `workdir` variable   |
                 |  in the "/var/www/../rna-seqlyze.wsgi" script must both    |
                 |  be set to the directory specified here!                   |
                 + ---------------------------------------------------------- +

Documentation:

    The 'workdir' holds

        - configuration files           \*.ini
        - the application database      rnaseqlyze.db
        - log files                     \*.log
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
    database.  If the command creates the sqlite database file, it changes it's
    unix access mode to (octal) 0664 and the group membership is changed to
    <group>.  <group> can be confgured in 'rnaseqlyze.ini'.  If the command
    creates the <workdir>, it changes it's unix access mode to (octal) 0775 and
    the group membership is also changed to <group>.  The command changes the
    unix access mode and group membership of all .log files inside the workdir
    to (octal) 0664 and <group>.
"""

import logging
log = logging.getLogger(__name__)

import os, sys, grp, shutil

import pkg_resources
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

import rnaseqlyze
import rnaseqlyze.web
import rnaseqlyze.worker
from rnaseqlyze.core.orm import Entity

Session = sessionmaker()

def main():

    import docopt
    opts = docopt.docopt(__doc__)

    workdir = os.path.abspath(opts['<workdir>'])

    # create the workdir if it does not exist
    wd_created = False
    if not os.path.isdir(workdir):
        if os.path.exists(workdir):
            log.error("not a directory: '%s'" % workdir)
            sys.exit(1)
        log.info("creating workdir '%s'" % workdir)
        os.makedirs(workdir)
        wd_created = True

    # create each config file that does not exist
    for pkg in rnaseqlyze, rnaseqlyze.web, rnaseqlyze.worker:

        # determine the destination file name
        conf_name = pkg.__name__.split('.')[-1] + ".ini"
        conf_path = os.path.join(workdir, conf_name)
        if os.path.exists(conf_path):
            continue

        # determine the source file name
        if pkg.__name__ == "rnaseqlyze":
            # for the core package there is currently
            # only one config file template
            ini = 'rnaseqlyze.ini'
        else:
            # for the non-core packages, take the desired version
            ini = opts['--development'] and "development.ini" or "production.ini"
        # get the file as a resource stream, which works even
        # if the distribution if installed as a zipped .egg
        req = pkg_resources.Requirement.parse(pkg.project_name)
        res = pkg_resources.resource_stream(req, ini)
        log.info("creating config file '%s'" % conf_name)
        shutil.copyfileobj(res, open(conf_path, "w"))

    # init rnaseqlyze configuration -- creates all .log files
    rnaseqlyze.configure(workdir)

    # set proper permissions on the log files
    for name in os.listdir(workdir):
        if name.endswith('.log'):
            path = os.path.join(workdir, name)
            log.info("adjusting permissions on '%s'" % name)
            os.chmod(path, 0664)
            os.chown(path, -1, grp.getgrnam(rnaseqlyze.group).gr_gid)

    # delayed because 'rnaseqlyze.group' was
    # not known before calling rnaseqlyze.configure() above
    if wd_created:
        log.info("adjusting permissions on '%s'" % workdir)
        # change permission bits
        os.chmod(workdir, 0775)
        # change group membership
        os.chown(workdir, -1, grp.getgrnam(rnaseqlyze.group).gr_gid)

    # get the database file path
    db_path = rnaseqlyze.db_url.split(":", 1)[1]

    # remove the databse file
    # if it exists and --recreatedb is given
    if os.path.exists(db_path) and opts['--recreatedb']:
        log.info("removing existing database file '%s'" %
                                                   db_path.split('/')[-1])
        os.unlink(db_path)

    # create the database if it doesn't exist
    if not os.path.exists(db_path):

        log.info("recreating database '%s'" % rnaseqlyze.db_url)

        # create sqlalchemy db engine
        engine = create_engine(rnaseqlyze.db_url)

        # create the file and initialize the schema
        with engine.begin() as conn:
            Entity.metadata.create_all(conn)

        log.info("adjusting permissions on database file")

        # change permission bits
        os.chmod(db_path, 0664)
        # change group membership
        os.chown(db_path, -1, grp.getgrnam(rnaseqlyze.group).gr_gid)

        log.info("initializing organism cache")

        # initialize UCSC Browser list of organisms
        from rnaseqlyze import org_cache
        with engine.begin() as conn:
            session = Session(bind=conn)
            org_cache.refresh(session)
            session.commit()

    log.info("workdir initialized")
