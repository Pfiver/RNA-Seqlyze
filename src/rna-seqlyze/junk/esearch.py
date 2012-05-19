#!/usr/bin/python
# encoding: utf-8

from Bio import Entrez
Entrez.email = "patrick.pfeifer@students.fhnw.ch"

import sys
handle = Entrez.esearch(*sys.argv[1:])
record = Entrez.read(handle)
print record
