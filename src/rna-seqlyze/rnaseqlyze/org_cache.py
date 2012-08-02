"""
RNA-Seqlyze keeps a cache of organisms available
in the UCSC Browser. In addition to that, for each of
those organisms, the matching refseq accession is cached.
"""
import logging
log = logging.getLogger(__name__)

import csv
import difflib
from pkg_resources import resource_stream

from rnaseqlyze import ucscbrowser
from rnaseqlyze.core.entities import UCSCOrganism

prokaryotes_tsv = "refseq-data/prokaryotes.txt"

def refresh(db_session):
    """
    Refresh the organism cache.

    The cache is initialized from the list of organisms available in the
    UCSC genome browser. A list of rnaseqlyze.orm.UCSCOrganism's is
    retrieved by calling rnaseqlyze.ucscbrowser.get_org_list().

    The retrieved list is not ready to be .add()ed to the :param:db_session
    however, because the objects' primary keys, the refseq accession,
    are still missing.

    Those are determined by parsing the list of gomplete genomes available in
    the ncbi "genome" database, which is stored in

        rnaseqlyze/refseq-data/prokaryotes.txt

    The file was retrieved from

        ftp://ftp.ncbi.nih.gov/genomes/GENOME_REPORTS/prokaryotes.txt

    on Mon, 02 Jul 2012.

    Once found, the rnaseqlyze.orm.UCSCOrganism objects are updated with the
    refseq accessions and .add()ed to the passed :param:db_session.
    """

    organisms = ucscbrowser.get_org_list()
    accessions = get_accessions()

    for org in organisms:
        ot = org.title
        for gt, acc in accessions:
            if ot == gt:
                org.acc = acc
                db_session.add(org)

    not_found = set(organisms) - \
                set(db_session.query(UCSCOrganism).all())

    # try fuzzy-matching the title
    # of those organisms that were not found
    for org in not_found:
        ot = org.title
        best_ratio = 0
        best_match = None
        for gt, acc in accessions:
            ratio = difflib.SequenceMatcher(None, ot, gt).ratio()
            if ratio > best_ratio:
                best_ratio = ratio
                best_match = acc
                best_match_t = gt

        if best_ratio > 0.8:
            if db_session.query(UCSCOrganism).get(best_match):
                log.debug(("NOT using '{match}' for '{org}'"
                           " despite match ratio of: {ratio}").format(
                               match=best_match_t, org=ot, ratio=best_ratio))
            else:
                log.info(("using '{match}' for '{org}'"
                          " match ratio: {ratio}").format(
                               match=best_match_t, org=ot, ratio=best_ratio))
                org.acc = best_match
                db_session.add(org)
        else:
            log.warn(("'{org}' not found in NCBI 'genome' database"
                      " (best match ratio only {ratio})").format(
                          org=ot, ratio=best_ratio))

    # make sure that that regular expression in views.post() that translates
    # the 'org_accession' from 'title (db/accession)', as generated in
    # rnaseqlyze.create.js, back to 'accession' doesn't fail
    for org in db_session.query(UCSCOrganism).all():
        if any(needle in heystack
                    for needle in '()'
                    for heystack in (org.db, org.acc, org.title)):
            log.warn("Droping organism with parentesis"
                     " to avoid problems in parsing auto"
                     "completed form input in views.post()")
            db_session.expunge(org)

def get_accessions():

    data_file = resource_stream(__name__, prokaryotes_tsv)
    reader = csv.reader(data_file, delimiter='\t')
    headings = reader.next()
    colnums = dict(zip(headings, map(headings.index, headings)))
    # -> { '#Organism/Name': 0, ..., 'Chromosomes/RefSeq': 7 }
    ret = []
    for cols in reader:
        if cols[colnums['Chromosomes/RefSeq']] == '-':
            continue
        ret.append((cols[colnums['#Organism/Name']],
                    cols[colnums['Chromosomes/RefSeq']]))
    return ret
