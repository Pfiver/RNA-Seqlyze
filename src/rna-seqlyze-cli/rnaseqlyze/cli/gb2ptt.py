import sys

from rnaseqlyze.gb2ptt import gb2ptt

def main(argv=sys.argv):
    gb2ptt(argv[1], sys.stdout)
