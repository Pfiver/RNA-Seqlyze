#!/usr/bin/python

########################################################################
#
# GEO query
#

GSE_acc = "GSE18630"

from string import Template
GSE_MINiML = Template("http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=$GSE_acc&form=xml&targ=self&view=full")

import urllib
from lxml import etree
#tree = etree.parse(urllib.urlopen(GSE_MINiML.substitute(GSE_acc=GSE_acc)))
tree = etree.parse(open(GSE_acc+".xml"))
find_samples = etree.ETXPath("/{%(ns)s}MINiML/{%(ns)s}Series[@iid=\"%(GSE)s\"]/{%(ns)s}Relation[@type=\"SRA\"]" % {
	"ns": tree.getroot().nsmap[None],
	"GSE": GSE_acc } )

SRPs = find_samples(tree)

if len(SRPs) != 1:
	raise Exception("unexpected number of SRA relations: %d", len(SRPs))

SRP_url = SRPs[0].attrib["target"]

import urlparse
elems = urlparse.urlparse(SRP_url)
query = urlparse.parse_qs(elems.query)

SRP_acc = query["term"][0]
print "SRA Platform link:", SRP_acc

########################################################################
#
# SRA query (SRP -> SRRs)
#

from Bio import Entrez
Entrez.email = "patrick.pfeifer@students.fhnw.ch"

handle = Entrez.esearch(db="sra", term=SRP_acc + "[Accession]")
record = Entrez.read(handle)

IDs = record["IdList"]

if len(IDs) < 1:
	raise Exception("unexpected number of SRA IDs: %d", len(IDs))

for ID in IDs:
	handle = Entrez.efetch(db="sra", rettype="xml", id=ID)
	open("%s-%s.xml" % (SRP_acc, ID), "w").writelines(handle)



# EXPERIMENT_PACKAGE_SET EXPERIMENT_PACKAGE RUN_SET RUN @accession
