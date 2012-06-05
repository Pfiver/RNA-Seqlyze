# encoding: utf-8

import sys

from rnaseqlyze import galaxy

def main(argv=sys.argv):
    galaxy.upload(argv[1])
