# encoding: utf-8

import sys

from rnaseqlyze.transterm import run

def main(argv=sys.argv):
    run(*argv[1:])
