"""\
RNA-Seqlyze DB-Init

Delete (!) & re-initialize the global database.

The command takes no arguments. The database is expected to be an sqlite
database and the path to the database file should be set in 'rnaseqlyze.ini'.
Afte creating the file, it's unix access mode is changed to (octal) 0664, which
means user/group=read/write,others=read and the group membership is changed to
<web_group>. <web_group> should be confgured in 'rnaseqlyze.ini' as well.
"""

import os, sys, grp

from sqlalchemy import create_engine

import rnaseqlyze
from rnaseqlyze.core.orm import Entity

def main():

    if len(sys.argv) > 1:
        print __doc__
        return

    # get the file path
    db_path = rnaseqlyze.db_url.split(":", 1)[1]
    # remove the file
    os.unlink(db_path)
    # recreate it & create all necessary tables
    Entity.metadata.create_all(create_engine(rnaseqlyze.db_url))
    # change parmission to user/group=read/write,others=read
    os.chmod(db_path, 0664)
    # change group ownership to web_server group
    os.chown(db_path, -1, grp.getgrnam(rnaseqlyze.web_group).gr_gid)

    print "initialized database:", db_path
