"""
RNA-Seqlyze keeps a cache of organisms available
in the UCSC Browser. In addition to that, for each of
those organisms, the matching "nucleaotide" database gi's,
accession numbers ('Caption') and 'Title's are cached.
"""

from Bio import Entrez

from rnaseqlyze import ucscbrowser

from rnaseqlyze.core.orm import NucleotideSummary

def refresh(db_session):
    """
    Refresh the organism cache.

    The cache is initialized from the list of organisms available in the
    UCSC genome browser. The list of rnaseqlyze.orm.UCSCOrganism's is
    retrieved by calling rnaseqlyze.ucscbrowser.get_org_list().
    This list is .add()ed to the passed :param:db_session.

    Then, for each organism, the following "esearch.fcgi" query is then made
    on the NCBI "nucleotide" database:

        db=nucleotide&term="<UCSC Organism Title>"[Title] AND refseq[Filter]

    The returned gi's are collected and then the "DocSummaries" are retrieved
    using an "esummary.fcgi" query.

    Finally, the 'Gi', 'Title' and 'Caption' fields of these "DocSummaries"
    are saved in the rnaseqlyze.orm.NucleaotideSummary table.
    """

    organisms = ucscbrowser.get_org_list()

    # fetch corresponding 'gi' id's
    # in batches of 25, in order not to exceed the maximum url length
    gis = []
    batch_size = 25
    for i in range(0, len(organisms), batch_size):
        titles = ' OR '.join('"%s"[Title]' % org.title
                             for org in organisms[i:i+batch_size])
        handle = Entrez.esearch(db="nucleotide", retmax=1<<30,
                                term='(%s) AND refseq[Filter]' % titles)
        result = Entrez.read(handle)
        gis.extend(result["IdList"])

    handle = Entrez.epost(db="nucleotide", id=",".join(gis))
    result = Entrez.read(handle)

    webEnv = result["WebEnv"]
    queryKey = result["QueryKey"]
    handle = Entrez.esummary(db="nucleotide", webenv=webEnv, query_key=queryKey)

    result = Entrez.read(handle)

    docsums = []
    for sum in result:
        docsums.append(NucleotideSummary(
            gi=sum['Gi'], title=sum['Title'], caption=sum['Caption']))

    db_session.add_all(organisms)
    db_session.add_all(docsums)
