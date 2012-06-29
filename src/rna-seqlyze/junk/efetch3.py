#!/usr/bin/python

def main():

    db = "genome"
    rettype = "text"

    from Bio import Entrez
    Entrez.email = "patrick.pfeifer@students.fhnw.ch"

    import sys
    args = sys.argv[1:]

    if len(args) == 0:

        import json
        import urllib2
    #    genomes = json.load(urllib2.urlopen(
    #        "http://archaea.ucsc.edu/wp-content/data/archaealViruses_genomes.json"))

        import sys
        genomes = json.load(sys.stdin)

        for genome in genomes:
            for child in genome['children']:
                for grandchild in child['children']:
                    print grandchild['data']['title']
    else:

        from pprint import pprint

        handle = Entrez.esummary(db="nucleotide", id=",".join(args))
        pprint(Entrez.read(handle))
        return

        if False:
            ids=[]
            for title in args:
                import sys
                handle = Entrez.esearch(db=db,
                        term=title + "[Organism]")
                record = Entrez.read(handle)
                ids.append((title, record["IdList"]
                                   and record["IdList"][0] or None))
            print ids
        else:
            handle = Entrez.esearch(db=db,
                        retmax=1<<30, term=" OR ".join(
                            title + "[Organism]" for title in args))
            record = Entrez.read(handle)
            print record["IdList"]
    return
main()
