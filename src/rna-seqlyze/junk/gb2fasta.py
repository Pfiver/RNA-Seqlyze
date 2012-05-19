#!/usr/bin/python

import sys
import Bio.SeqIO
Bio.SeqIO.write(
        Bio.SeqIO.parse(
                open(sys.argv[1]),
                "genbank"),
        sys.stdout,
        "fasta")