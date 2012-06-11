from __future__ import print_function

import sys
from sqlalchemy import create_engine

import rnaseqlyze
from rnaseqlyze.cli import DBSession
from rnaseqlyze.core.orm import Entity

def main(argv=sys.argv):
    if len(argv) > 1:
        print("%s: ignoring arguments" % argv[0], file=sys.stderr)
    engine = create_engine(rnaseqlyze.db_url)
    DBSession.configure(bind=engine)
    Entity.metadata.create_all(engine)
