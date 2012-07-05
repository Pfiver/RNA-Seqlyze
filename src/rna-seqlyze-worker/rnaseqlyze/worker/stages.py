"""
RNA-Seqlyze Worker Stages

    -- **this** is where things are actually getting done! :-)
"""

import os
from os.path import join, exists, isdir
from subprocess import Popen, PIPE
from StringIO import StringIO
from urllib import quote

from Bio import SeqIO
from Bio.SeqFeature import \
        SeqFeature, FeatureLocation, ExactPosition

from rnaseqlyze import efetch
from rnaseqlyze import galaxy
from rnaseqlyze import ucscbrowser
from rnaseqlyze import transterm, gb2ptt
from rnaseqlyze.ucscbrowser import BAMTrack, BigWigTrack, BigBedTrack
from rnaseqlyze.core.entities import GalaxyDataset

_stages = []
def stage(method):
    """
    Just a small helper to collect the stages in the order defined.

    To add a new stage, simply add a method to :class:`~WorkerStages`.
    It will be automatically executed for all new analyses.
    """
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

    ######################################################################
    # Utility Methods & Properties

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

    @property
    def coverage_name(self):
        return self.bam_name + "-coverage"

    @property
    def hpterms_name(self):
        return self.analysis.org_accession + "-hairpin-terminators"

    ##################################################################
    # Stages

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
    def fetch_srr(self):
        # download if not already in cache
        if not os.path.exists(self.analysis.rnaseq_run.sra_path):
            self.log.debug("transfering input file from sra")
            self.analysis.rnaseq_run.download()
            self.log.debug("done")

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
    def bam_to_wiggle_to_bigwig(self):
        if exists(join(self.data_dir, "tophat-output", "accepted_hits.bam")) \
           and not exists(join(self.data_dir, "coverage.bigwig")):
            os.chdir(self.data_dir)
            self.log_cmd("bam_to_wiggle.py", "-o", "coverage.bigwig",
                                        "tophat-output/accepted_hits.bam")

    @stage
    def transterm_hp(self):

        acc = self.analysis.org_accession
        ptt_path = join(self.gb_data_dir, acc + ".ptt")
        dummy_ptt_path = join(self.gb_data_dir, "chr.ptt")

        if not exists(dummy_ptt_path):
            os.symlink(acc + ".ptt", dummy_ptt_path)

        if not exists(join(self.data_dir, ptt_path)):
            self.log.debug("converting %s to ptt" % 
                      self.analysis.genbankfile_path.split('/')[-1])
            ptt_file = open(ptt_path, "w")
            gb_file = open(self.analysis.genbankfile_path)
            gb2ptt.gb2ptt(gb_file, ptt_file)
            ptt_file.close()
            gb_file.close()

        if not exists(join(self.data_dir, "hpterminators.bed")):
            os.chdir(self.data_dir)
            self.log.debug("running transterm")

            tt_out = open("transterm_hp.out", "w+")
            transterm.run((self.analysis.fa_path,
                           dummy_ptt_path), out=tt_out, err=self.logfile)
            tt_out.seek(0)
            bed_file = open("hpterminators.bed", "w")
            transterm.fa2bed(tt_out, bed_file)
            bed_file.close()
            tt_out.close()

            record = SeqIO.parse(open(self.analysis.fa_path),"fasta").next()
            self.log.debug("'chromosome' length: %d", len(record.seq))

            self.log.debug("running bedToBigBed")

            chrs = open("chrom.sizes", "w")
            chrs.write("chr %d" % len(record.seq))
            chrs.close()
            self.log_cmd("bedToBigBed",
                    "hpterminators.bed", "chrom.sizes", "hpterminators.bigbed")


    @stage
    def galaxy_upload(self):

        # FIXME: names are not unique on galaxy:
        # is "%s_%s" % (srr_name, self.analysis.org_accession) good enough ?

        if not self.analysis.galaxy_bam:
            bam_path = join(self.data_dir, 
                             "tophat-output", "accepted_hits.bam")
            self.log.info("uploading accepted_hits.bam to galaxy")
            self.analysis.galaxy_bam = GalaxyDataset(
                    galaxy.upload(open(bam_path), self.bam_name))
            self.log.info("...done - id: %s" % self.analysis.galaxy_bam.id)
            self.session.commit()
            return

        if not self.analysis.galaxy_coverage:
            coverage_path = join(self.data_dir, "coverage.bigwig")
            self.log.info("uploading coverage.bigwig to galaxy")
            self.analysis.galaxy_coverage = GalaxyDataset(
                    galaxy.upload(open(coverage_path), self.coverage_name))
            self.log.info("...done - id: %s" % self.analysis.galaxy_coverage.id)
            self.session.commit()

        if not self.analysis.galaxy_hpterms:
            hpterms_path = join(self.data_dir, "hpterminators.bigbed")
            self.log.info("uploading hpterminators.bigbed galaxy")
            self.analysis.galaxy_hpterms = GalaxyDataset(
                    galaxy.upload(open(hpterms_path), self.hpterms_name))
            self.log.info("...done - id: %s" % self.analysis.galaxy_hpterms.id)
            self.session.commit()


    @stage
    def create_and_upload_htg_custom_text(self):
        """
        ``hgt.customText`` is a paremeter of the UCSC
        "hgTracks" genome browser application that makes it
        possible to share "cutom tracks" via a url.

        The value of the ``hgt.customText`` parameter is itself
        an URL. The shareable "custom tracks url" is therefore an
        URL that containes another URL. The other url must be "escaped"
        for this to work. That actually happens in
        :meth:`~rnaseqlyze.core.analysis.AnalysisMixins.hg_url`.
        
        The details are explained here:
        http://genome.ucsc.edu/goldenPath/help/customTrack.html#SHARE
        """

        if self.analysis.galaxy_hgtext:
            return

        tracks = []

        # bam track (mapping)
        bam_url = "https://" + galaxy.hostname \
                    + galaxy.ucsc_bam_path_template \
                        .format(dataset=self.analysis.galaxy_bam.id)
        tracks.append(BAMTrack(url=bam_url,
                               name="RNA-Seqlyze | %s" % self.bam_name))

        # bigwig track (coverage)
        coverage_url = "https://" + galaxy.hostname \
                        + galaxy.dataset_display_url_template \
                            .format(dataset=self.analysis.galaxy_coverage.id)
        tracks.append(BigWigTrack(url=coverage_url,
                                  name="RNA-Seqlyze | %s" % self.coverage_name))

        # bigbed track (terminators)
        hpterms_url = "https://" + galaxy.hostname \
                        + galaxy.dataset_display_url_template \
                            .format(dataset=self.analysis.galaxy_hpterms.id)
        tracks.append(BigBedTrack(url=hpterms_url,
                                  name="RNA-Seqlyze | %s" % self.hpterms_name))

        track_file = StringIO()
        track_file.write('\n'.join(tracks))
        track_file.seek(0)
        self.analysis.galaxy_hg_custom_id = \
                    galaxy.upload(track_file, self.bam_name + ".txt")
        self.session.commit()

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
del _stages

class UnknownInputfileTypeException(Exception):
    pass
