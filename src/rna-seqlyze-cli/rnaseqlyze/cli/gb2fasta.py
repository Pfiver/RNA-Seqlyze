import sys

import Bio.SeqIO

def main(argv=sys.argv):
    Bio.SeqIO.write(
        Bio.SeqIO.parse(
            open(argv[1]),
            "genbank"),
        sys.stdout,
        "fasta")
