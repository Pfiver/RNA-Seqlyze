"""\
RNA-Seqlyze gb2ptt

Convert a genbank file to ptt (protein table) format

Usage:
    rnas-gb2ptt <input.gb> <output.ptt>

If <input.gb> is '-', use 'sys.stdin, if <output.ptt> is '-', use 'sys.stdout'.
"""
import sys, logging
from rnaseqlyze.gb2ptt import gb2ptt

def main():

    if len(sys.argv) < 2 or sys.argv[1] in ('-h', '--help'):
        print __doc__
        return

    inputfile = sys.argv[1] == '-' and sys.stdin or open(sys.argv[1])
    outputfile = sys.argv[2] == '-' and sys.stdout or open(sys.argv[2], "w")

    loggin.basicConfig(level=logging.NOTSET) # logs to stderr

    gb2ptt(inputfile, outputfile)
