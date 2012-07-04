import os
from os.path import join, exists, isdir
from subprocess import Popen, PIPE
from StringIO import StringIO
from urllib import quote

from Bio import SeqIO
from Bio.SeqFeature import \
        SeqFeature, FeatureLocation, ExactPosition

from rnaseqlyze import galaxy
from rnaseqlyze import ucscbrowser

_stages = []
def stage(method):
    _stages.append(method)
    return method

class WorkerStages(object):

    """
    Available attributes:

        - self.analysis
        - self.data_dir = self.analysis.data_dir
        - self.gb_data_dir = self.analysis.gb_data_dir
        - self.input_data_dir = self.analysis.input_data_dir
        - self.session = DBSession(bind=create_engine(rnaseqlyze.db_url))

    .. note::
        After changing one of the attributes of the self.analysis object,
        **allways** **immediately** call self.session.commit(). Otherwise
        the database will stay locked and the web frontend can't update the ui.

    """

    def log_cmd(self, *cmd):
        proc = Popen(cmd, stdout=self.logfile, stderr=self.logfile)
        proc.wait()
        if proc.returncode != 0:
            raise Exception("%s failed" % (cmd,))

    @property
    def srr_name(self):
        if self.analysis.inputfile_name: # uploaded input file
            return self.analysis.inputfile_name.rsplit(".", 1)[0]
        else: # input is an SRR Indentifier
            return self.analysis.rnaseq_run.srr

    @property
    def bam_name(self):
        return "%s_%s" % (self.srr_name, self.analysis.org_accession)

    @stage
    def determine_inputfile_type(self):
        def getit(header):
            if header[0] == '@': return 'fastq'
            elif header[:8] == 'NCBI.sra': return 'sra'
            else: raise UnknownInputfileTypeException()
        self.analysis.inputfile_type = \
                getit(open(self.analysis.inputfile_path).read(8))
        self.session.commit()

    @stage
    def convert_input_file(self):
        fq_path = self.analysis.inputfile_fqpath
        if not exists(fq_path):
            os.chdir(self.input_data_dir)
            if self.analysis.inputfile_name:
                sra_name = self.analysis.inputfile_name
            else:
                sra_name = self.analysis.rnaseq_run.sra_name
            self.log.info("converting %s" % sra_name)
            self.log_cmd("fastq-dump", sra_name)

    @stage
    def fetch_gb(self):
        acc = self.analysis.org_accession
        self.analysis.create_gb_data_dir()
        out_path = join(self.gb_data_dir, acc + ".gb")
        if not exists(out_path):
            self.log.info("Fetching '%s' from entrez..." % acc)
            # TODO: do this earlier
            gb_id = efetch.get_nc_id(acc)
            efetch.fetch_nc_gb(gb_id, open(out_path, "w"))
            self.log.info("...done")

    @stage
    def gb2fasta(self):
        acc = self.analysis.org_accession
        gb_path = join(self.gb_data_dir, acc + ".gb")
        fa_path = join(self.gb_data_dir, acc + ".fa")
        if not exists(fa_path):
            self.log.info("Converting '%s' to fasta format..." % acc)
            record = SeqIO.parse(open(gb_path), "genbank").next()
            record.id = "chr" # required (!) for ucsc browser
            SeqIO.write(record, open(fa_path, "w"), "fasta")

    @stage
    def bowtie_build(self):
        acc = self.analysis.org_accession
        bt2_path = join(self.gb_data_dir, acc + ".1.bt2")
        if not exists(bt2_path):
            os.chdir(self.gb_data_dir)
            self.log.info("bowtie2-build %s" % acc)
            self.log_cmd("bowtie2-build", acc + ".fa", acc)

    @stage
    def tophat(self):
        acc = self.analysis.org_accession
        if not exists(join(self.data_dir,
                "tophat-output", "accepted_hits.bam")):
            os.chdir(self.data_dir)
            n_cpus = os.sysconf("SC_NPROCESSORS_ONLN")
            acc_path = join(self.gb_data_dir, acc)
            self.log_cmd("tophat", "-p", str(n_cpus),
                    "-o", "tophat-output",
                    acc_path, self.analysis.inputfile_fqpath)

    @stage
    def galaxy_upload(self):
        if self.analysis.galaxy_bam_id:
            return
        bam_path = join(self.data_dir, 
                         "tophat-output", "accepted_hits.bam")
        # FIXME: names are not unique on galaxy:
        # is "%s_%s" % (srr_name, self.analysis.org_accession) good enough ?
        bam_file = open(bam_path)
        self.log.info("uploading %s to galaxy server %s ..." % (
                                 bam_path,           galaxy.hostname))
        self.analysis.galaxy_bam_id = \
                galaxy.upload(bam_file, self.bam_name)
        self.log.info("done.")
        bam_file.close()
        self.session.commit()

    @stage
    def create_and_upload_galaxy_ucsc_bam_track(self):
        if self.analysis.galaxy_ucsc_bam_track_id:
            return
        bam_url = "https://" + galaxy.hostname \
                    + galaxy.ucsc_bam_path_template \
                        .format(dataset=self.analysis.galaxy_bam_id)
        track_line = ucscbrowser.bam_track_line_template \
                        .format(track_name="RNA-Seqlyze | %s" %
                                self.bam_name, big_data_url=bam_url)
        track_file = StringIO()
        track_file.write(track_line)
        track_file.seek(0)
        self.analysis.galaxy_ucsc_bam_track_id = \
                    galaxy.upload(track_file, self.bam_name + ".txt")

    @stage
    def create_genbank_file(self):
        """
        Greate a genbank file containing

        For more documentation on how to create new features, visit

         - http://biopython.org/\\
                 DIST/docs/api/Bio.SeqRecord.SeqRecord-class.html#__getitem__
         - http://biopython.org/\\
                 DIST/docs/api/Bio.SeqFeature.SeqFeature-class.html

         - http://www.ebi.ac.uk/\\
                 embl/Documentation/FT_definitions/feature_table.html
        """

        class Operon(object):
            def __init__(self, **kwargs):
                self.__dict__.update(kwargs)
        
        operons = Operon(beg=0, end=100, strand=1),

        record = SeqIO.parse(open(
            self.analysis.genbankfile_path), "genbank").next()

        for i, oper in enumerate(operons):
            location = FeatureLocation(ExactPosition(oper.beg),
                                       ExactPosition(oper.end))
            record.features.append(
                SeqFeature(location,
                    type='mRNA',
                    strand=oper.strand,
                    qualifiers=dict(
                        note='putative',
                        operon='operon%d' % i)))

        record.features.sort(key=lambda f: f.location.start.position)

        xgb_file = open(self.analysis.xgenbankfile_path, "w")

        SeqIO.write(record, xgb_file, "genbank")

WorkerStages.members = _stages
del stage, _stages

class UnknownInputfileTypeException(Exception):
    pass
