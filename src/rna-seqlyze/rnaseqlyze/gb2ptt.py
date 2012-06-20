import sys

from Bio import SeqIO
from Bio.SeqFeature import ExactPosition

debug = False

def gb2ptt(in_file, out_file):
    seq = SeqIO.parse(in_file, "genbank").next()

    # input:
        #      CDS             249..857
        #                      /locus_tag="SSO0001"
        #                      /note="Predicted membrane protein, conserved in archaea"
        #                      /codon_start=1
        #                      /transl_table=11
        #                      /product="hypothetical protein"
        #                      /protein_id="NP_341578.1"
        #                      /db_xref="GI:15896973"
        #                      /db_xref="GeneID:1455258"
        #                      /translation="MITEFLLKKKLEEHLSHVKEENTIYVTDLVRCPRRVRYESEYKE
        #                      LAISQVYAPSAILGDILHLGLESVLKGNFNAETEVETLREINVGGKVYKIKGRADAII
        #                      RNDNGKSIVIEIKTSRSDKGLPLIHHKMQLQIYLWLFSAEKGILVYITPDRIAEYEIN
        #                      EPLDEATIVRLAEDTIMLQNSPRFNWECKYCIFSVICPAKLT"

    # output:
        # Sulfolobus solfataricus P2 chromosome, complete genome - 1..2992245
        # 2978 proteins
        # Location        Strand  Length  PID     Gene    Synonym Code    COG     Product
        # 249..857        +       202     15896973        -       SSO0001 -       COG1468L        hypothetical protein

    import csv
    writer = csv.writer(out_file, delimiter='\t', lineterminator='\n')
    writer.writerow((seq.description,))
    writer.writerow(())
    writer.writerow(('Location', 'Strand', 'Length', 'PID', 'Gene', 'Synonym', 'Code', 'COG', 'Product'))

    import logging
    log = logging.getLogger(__name__)
    log.addHandler(logging.StreamHandler(sys.stderr))
    log.setLevel([logging.INFO, logging.DEBUG][debug])

    n=0
    for f in seq.features:

        n+=1
        if debug and n > 10:
            break

        if f.type != 'CDS':
            continue

        if type(f.location.start) != ExactPosition \
           or type(f.location.end) != ExactPosition:
            log.info("skipping non-exact location '%s' in '%s'" % (f.location, f.type))
            continue

        _len = f.location.end.position - f.location.start.position
        if _len < 0:
            _len = len(seq.seq) - f.location.start.position + f.location.end.position
        if _len % 3:
            log.info("implausible feature length (%d) in '%s'" % (_len, f.type))
        _len //= 3 # integer division
        _len -= 1 # omit stop codon

        xrefs = dict(map(lambda s: s.split(':'), f.qualifiers['db_xref']))
        for r in 'GI', 'GeneID':
            if r not in xrefs:
                xrefs[r] = '-'

        for q in 'gene', 'product', 'locus_tag':
            if q not in f.qualifiers:
                f.qualifiers[q] = '-'

        writer.writerow((
            # convert between biopython (0-based, incl:excl) and genbank (1-based, incl:incl) positions
            "%d..%d" % (f.location.start.position+1, f.location.end.position), ['-', '+'][(f.strand + 1) / 2], _len,
            xrefs['GI'], f.qualifiers['gene'][0], f.qualifiers['locus_tag'][0], '-', '-', f.qualifiers['product'][0]))

    log.info("wrote %d rows" % n)
