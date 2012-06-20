import sys

from rnaseqlyze import galaxy

def main(argv=sys.argv):
    print galaxy.upload(oen(argv[1]), os.path.basename(argv[1]))
