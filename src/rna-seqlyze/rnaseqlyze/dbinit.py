from __future__ import print_function

import sys
from sqlalchemy import create_engine

import rnaseqlyze
from rnaseqlyze.core.orm import Entity

def main(argv=sys.argv):

    if len(argv) > 1:
        print("%s: ignoring arguments" % argv[0], file=sys.stderr)

    Entity.metadata.create_all(create_engine(rnaseqlyze.db_url))
