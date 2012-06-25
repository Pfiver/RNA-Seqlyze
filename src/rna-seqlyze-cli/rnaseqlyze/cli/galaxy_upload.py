"""
Galaxy-Upload

Usage:
    galaxy-upload <local_file>
"""

import sys
import rnaseqlyze

def main(argv=sys.argv):
    print rnaseqlyze.galaxy.upload(
            oen(argv[1]), os.path.basename(argv[1]))
