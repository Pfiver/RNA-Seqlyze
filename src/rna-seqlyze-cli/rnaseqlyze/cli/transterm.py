"""\
RNA-Seqlyze transterm

Calls the transterm program
with an additional "-p <path to>/expterm.dat" argument.

Usage:
    rnas-transterm [--] <transterm arguments> ...
    rnas-transterm -h|--help
"""
import sys
from rnaseqlyze.transterm import run

def main():

    if len(sys.argv) > 1 and sys.argv[1] == '--':
        sys.argv.pop(1)
    elif len(sys.argv) < 2 or sys.argv[1] in ('-h', '--help'):
        print __doc__
        return

    run(sys.argv[1:])
