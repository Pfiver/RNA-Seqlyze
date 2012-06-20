import sys, logging

def main(argv=sys.argv):

    loggin.basicConfig(level=logging.NOTSET)

    from rnaseqlyze.gb2ptt import gb2ptt
    gb2ptt(open(argv[1]), sys.stdout)
