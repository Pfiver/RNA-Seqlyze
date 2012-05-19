#!/usr/bin/python

#db = "nuccore"
db = "gds"
#rettype = "gb"
rettype = "xml"
#accession = "NC_002754"
#accession = "SRP001461"
accession = "GSE18630"

from Bio import Entrez
Entrez.email = "patrick.pfeifer@students.fhnw.ch"

handle = Entrez.esearch(db=db, term=accession + "[Accession]")
record = Entrez.read(handle)

IDs = record["IdList"]

if len(IDs) != 1:
	raise Exception("unexpected reply from Entrez: got %d IDs: %s" % (len(IDs),IDs))

handle = Entrez.efetch(db=db, id=IDs[0], rettype=rettype, retmode="text")

open(accession + "." + rettype, "w").writelines(handle)
