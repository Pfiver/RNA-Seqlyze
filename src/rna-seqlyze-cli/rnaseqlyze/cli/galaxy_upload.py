"""\
RNA-Seqlyze Galaxy-Upload

Usage:
    rnas-galaxy-upload <local_file>
"""
import sys, os
from rnaseqlyze import galaxy

def main():

    if len(sys.argv) < 2 or sys.argv[1] in ('-h', '--help'):
        print __doc__
        return

    print galaxy.upload(open(sys.argv[1]), os.path.basename(sys.argv[1]))
