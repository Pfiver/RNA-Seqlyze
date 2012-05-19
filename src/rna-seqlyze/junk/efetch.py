#!/usr/bin/python
# encoding: utf-8

from Bio import Entrez
Entrez.email = "patrick.pfeifer@students.fhnw.ch"

import sys
handle = Entrez.efetch(db=sys.argv[1], id=sys.argv[2], rettype="xml")
print handle.read()
