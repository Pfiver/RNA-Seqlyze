from __future__ import print_function

import sys
from sqlalchemy import create_engine

import rnaseqlyze
from rnaseqlyze.core.orm import Entity

def main(argv=sys.argv):
    import os, grp
    if len(argv) > 1:
        print("%s: ignoring arguments" % argv[0], file=sys.stderr)
    # get path
    db_path = rnaseqlyze.db_url.split(":", 1)[1]
    # remove
    os.unlink(db_path)
    #recreate
    Entity.metadata.create_all(create_engine(rnaseqlyze.db_url))
    # change parmission to user/group=read/write,others=read
    os.chmod(db_path, 0664)
    # change group ownership to web_server group
    os.chown(db_path, -1, grp.getgrnam(rnaseqlyze.web_group).gr_gid)
