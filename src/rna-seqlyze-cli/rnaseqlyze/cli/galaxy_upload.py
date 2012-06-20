import sys

from rnaseqlyze import galaxy

def main(argv=sys.argv):
    galaxy.upload(argv[1], os.path.basename(argv[1]))
