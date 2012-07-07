"""
RNA-Seqlyze Worker Stages

    -- **this** is where things are actually getting done! :-)
"""

import os
from os.path import join, exists, isdir, relpath
from subprocess import Popen, PIPE
from StringIO import StringIO
from urllib import quote

import pysam

from Bio import SeqIO
from Bio.SeqFeature import \
        SeqFeature, FeatureLocation, ExactPosition

from rnaseqlyze import efetch
from rnaseqlyze import galaxy
from rnaseqlyze import ucscbrowser
from rnaseqlyze import transterm, gb2ptt
from rnaseqlyze.ucscbrowser import BAMTrack, BigWigTrack, BigBedTrack
from rnaseqlyze.core.entities import GalaxyDataset


class Operon(object):
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

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
        - self.session = DBSession(bind=create_engine(rnaseqlyze.db_url))

    .. note::
        After changing one of the attributes of the self.analysis object,
        **allways** **immediately** call self.session.commit(). Otherwise
        the database will stay locked and the web frontend can't update the ui.

    """

    ############################################################################
    # Utility Methods & Properties

    def log_cmd(self, *cmd):
        proc = Popen(cmd, stdout=self.logfile, stderr=self.logfile)
        proc.wait()
        if proc.returncode != 0:
            raise Exception("%s failed" % (cmd,))

    @property
    def srr_name(self):
        return self.analysis.inputfile_base_name

    @property
    def bam_name(self):
        return "%s %s Mapping" % (self.genbank_record.id, self.srr_name)

    @property
    def coverage_name(self):
        return "%s %s Coverage" % (self.genbank_record.id, self.srr_name)

    @property
    def hp_terms_name(self):
        return "%s Hairpin Terminators" % (self.genbank_record.id,)

    @property
    def pr_operons_name(self):
        return "%s Predicted Operons" % (self.genbank_record.id,)

    #
    ############################################################################
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
        if not exists(self.analysis.inputfile_fq_path):
            self.log.debug("creating %s" % self.analysis.inputfile_fq_path)
            os.chdir(self.analysis.input_data_dir)
            self.log_cmd("fastq-dump", self.analysis.inputfile_name)

    @stage
    def fetch_genbank_file(self):
        if not exists(self.analysis.genbankfile_path):
            if not exists(self.analysis.genbank_data_dir):
                os.makedirs(self.analysis.genbank_data_dir)
            self.log.info("Fetching '%s' from entrez..." %
                                     self.analysis.org_accession)
            gb_id = efetch.get_nc_id(self.analysis.org_accession)
            efetch.fetch_nc_gb(gb_id, open(self.analysis.genbankfile_path, "w"))
            self.log.info("...done")

    @stage
    def read_genbank_file(self):
        self.genbank_record = SeqIO.parse(open(self.analysis \
                                    .genbankfile_path), "genbank").next()
        ngenes = sum(1 for f in self.genbank_record.features
                       if f.type == 'gene')
        self.log.info("genbank file lists %d genes" % ngenes)

    @stage
    def gb2fasta(self):
        if not exists(self.analysis.genbankfile_fa_path):
            self.log.info("Converting '%s' to fasta format..." %
                    self.analysis.genbankfile_path)
            saved_id = record.id
            record.id = "chr" # required (!) for ucsc browser
            SeqIO.write(record, open(
                self.analysis.genbankfile_fa_path, "w"), "fasta")
            record.id = saved_id

    @stage
    def bowtie_build(self):
        if not exists(join(self.analysis.genbank_data_dir,
                           self.analysis.genbankfile_base_name + ".1.bt2")):
            os.chdir(self.analysis.genbank_data_dir)
            self.log_cmd("bowtie2-build",
                            self.analysis.genbankfile_fa_name,
                            self.analysis.genbankfile_base_name)

    @stage
    def tophat(self):
        if not exists(join(self.analysis.data_dir,
                           "tophat-output", "accepted_hits.bam")):
            os.chdir(self.analysis.data_dir)
            n_cpus = os.sysconf("SC_NPROCESSORS_ONLN")
            self.log_cmd("tophat", "-p", str(n_cpus),
                    "-o", "tophat-output",
                    relpath(join(self.analysis.genbank_data_dir,
                                 self.analysis.genbankfile_base_name)),
                    relpath(self.analysis.inputfile_fq_path))

    @stage
    def bam_to_bigwig(self):
        if exists(join(self.analysis.data_dir,
                       "tophat-output", "accepted_hits.bam")) \
           and not exists(join(self.analysis.data_dir, "coverage.bigwig")):
            os.chdir(self.analysis.data_dir)
            # the script automatically converts it's
            # output to bigwig if it finds kent's wigToBigWig
            self.log_cmd("bam_to_wiggle.py", "-o", "coverage.bigwig",
                                        "tophat-output/accepted_hits.bam")

    @stage
    def transterm_hp(self):
        ptt_name = self.genbank_record.id + ".ptt"
        ptt_path = join(self.analysis.genbank_data_dir, ptt_name)
        dummy_ptt_path = join(self.analysis.genbank_data_dir, "chr.ptt")

        if not exists(dummy_ptt_path):
            os.symlink(ptt_name, dummy_ptt_path)

        if not exists(join(self.analysis.data_dir, ptt_path)):
            self.log.debug("converting %s to ptt" % 
                                       self.analysis.genbankfile_name)
            ptt_file = open(ptt_path, "w")
            gb_file = open(self.analysis.genbankfile_path)
            gb2ptt.gb2ptt(gb_file, ptt_file)
            ptt_file.close()
            gb_file.close()

        if not exists(join(self.analysis.data_dir, "hp_terminators.bigbed")):
            os.chdir(self.analysis.data_dir)
            self.log.debug("running transterm")

            tt_out = open("transterm_hp.out", "w+")
            # --min-conf=n n is the cut-off confidence value,
            #              between 0 and 100, the default is 76
            tt_args = ("--min-conf=47",
                       self.analysis.genbankfile_fa_path, dummy_ptt_path)
            transterm.run(tt_args, out=tt_out, err=self.logfile)
            tt_out.seek(0)

            self.hp_terminators = list(transterm.iterator(tt_out))
            tt_out.seek(0)
            self.log.info("found {0} possible hairpin terminators"
                                     .format(len(self.hp_terminators)))

            bed_file = open("hp_terminators.bed", "w")
            transterm.tt2bed(tt_out, bed_file)
            bed_file.close()
            tt_out.close()

            self.log.debug("running bedToBigBed")

            chrs = open("chrom.sizes", "w")
            chrs.write("chr %d" % len(self.genbank_record.seq))
            chrs.close()
            self.log_cmd("bedToBigBed", "hp_terminators.bed",
                             "chrom.sizes", "hp_terminators.bigbed")

    @stage
    def predict_operons(self):
        bam_path = join(join(self.analysis.data_dir,
                        "tophat-output", "accepted_hits.bam"))
        if not exists(bam_path + ".bai"):
            pysam.index(bam_file)

        self.max = 0
        self.covered = 0
        self.coverage = [0] * len(self.genbank_record.seq)

        sam_reader = pysam.Samfile(bam_path, "rb")
        chrom, length = sam_reader.references[0], sam_reader.lengths[0]

        assert chrom == "chr" and length == len(self.genbank_record.seq), (
                "Something went badly wrong"
                " -- the bam track or genbank file cold be corrupted...")

        for base in sam_reader.pileup(chrom, 0, length):
            self.covered += 1
            if base.n > self.max:
                self.max = base.n
            self.coverage[base.pos] = base.n

        if not self.covered:
            raise Exception("Not a valid bam file")

        self.log.debug("maximum coverage: %d" % self.max)
        self.log.debug("number of bases covered by short reads: %d/%d" % (
                                    self.covered, len(self.genbank_record.seq)))


        # available here:
        # ---------------
        #
        # - self.genbank_record: Biopython SeqIO.parse()d genbank file
        #
        # - self.coverage: [n,n,n,n,...] / len = len(self.genbank_record.seq)
        #
        # - self.hpterminators: ((id, begin, end, strand, confidence), ...)
        #                         str, str, str, str (1/-), int
        #

        self.operons = Operon(begin=0, end=100, strand=1, confidence=10),
        self.operons = Operon(begin=200, end=300, strand=1, confidence=50),
        self.operons = Operon(begin=400, end=500, strand=1, confidence=100),

    @stage
    def create_operon_track(self):
        track_name = "rna-seqlyze-operon_predictions"
        os.chdir(self.analysis.data_dir)
        if exists(track_name + ".bigbed"):
            return

        bed_file = open(track_name + ".bed", "w")
        for i, o in enumerate(self.operons):
            begin, end = str(o.begin), str(o.end)
            rgb_color = ','.join((str(100 - int(o.confidence)),)*3)
            print >> bed_file, '\t'.join((
                'chr', begin, end,
                'OPERON_%d' % i, str(o.confidence),
                '+' if o.strand > 0 else '-', begin, end, rgb_color
            ))
        bed_file.close()

        # chrom_sizes already generated during "transterm_hp"
        self.log_cmd("bedToBigBed", track_name + ".bed",
                         "chrom.sizes", track_name + ".bigbed")

    @stage
    def upload_track_data(self):

        # FIXME: names are not unique on galaxy:
        # is "%s_%s" % (srr_name, self.analysis.org_accession) good enough ?

        if not self.analysis.galaxy_bam:
            bam_path = join(self.analysis.data_dir, 
                             "tophat-output", "accepted_hits.bam")
            self.log.info("uploading accepted_hits.bam to galaxy")
            self.analysis.galaxy_bam = GalaxyDataset(
                id=galaxy.upload(open(bam_path), self.bam_name))
            self.log.info("...done - id: %s" % self.analysis.galaxy_bam.id)
            self.session.commit()

        if not self.analysis.galaxy_coverage:
            coverage_path = join(self.analysis.data_dir, "coverage.bigwig")
            self.log.info("uploading coverage.bigwig to galaxy")
            self.analysis.galaxy_coverage = GalaxyDataset(
                id=galaxy.upload(open(coverage_path), self.coverage_name))
            self.log.info("...done - id: %s" % self.analysis.galaxy_coverage.id)
            self.session.commit()

        if not self.analysis.galaxy_hp_terms:
            hp_terms_path = join(self.analysis.data_dir,
                                 "hp_terminators.bigbed")
            self.log.info("uploading hp_terminators.bigbed to galaxy")
            self.analysis.galaxy_hp_terms = GalaxyDataset(
                id=galaxy.upload(open(hp_terms_path), self.hp_terms_name))
            self.log.info("...done - id: %s" % self.analysis.galaxy_hp_terms.id)
            self.session.commit()

        if not self.analysis.galaxy_pr_operons:
            track_filename = "rna-seqlyze-operon_predictions.bigbed"
            pr_operons_path = join(self.analysis.data_dir, track_filename)
            self.log.info("uploading %s to galaxy" % track_filename)
            self.analysis.galaxy_pr_operons = GalaxyDataset(
                id=galaxy.upload(open(pr_operons_path), self.pr_operons_name))
            self.log.info("...done - id: %s" %
                                         self.analysis.galaxy_pr_operons.id)
            self.session.commit()

    @stage
    def create_and_upload_hg_text(self):
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

        if self.analysis.galaxy_hg_text:
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
        hp_terms_url = "https://" + galaxy.hostname \
                        + galaxy.dataset_display_url_template \
                            .format(dataset=self.analysis.galaxy_hp_terms.id)
        tracks.append(BigBedTrack(url=hp_terms_url,
                                  name="RNA-Seqlyze | %s" % self.hp_terms_name))

        # bigbed track (predicted operons)
        pr_operons_url = "https://" + galaxy.hostname \
                        + galaxy.dataset_display_url_template \
                            .format(dataset=self.analysis.galaxy_pr_operons.id)
        tracks.append(BigBedTrack(url=pr_operons_url,
                                  name="RNA-Seqlyze | %s" %
                                                      self.pr_operons_name))

        track_file = StringIO()
        track_file.write('\n'.join(tracks))
        track_file.seek(0)
        self.analysis.galaxy_hg_text = GalaxyDataset(
                    id=galaxy.upload(track_file,
                                     "UCSC Tracks Analysis%d.txt" %
                                                          self.analysis.id))
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
        if exists(self.analysis.xgenbankfile_path):
            return

        record = SeqIO.parse(open(
            self.analysis.genbankfile_path), "genbank").next()

        for i, o in enumerate(self.operons):
            location = FeatureLocation(ExactPosition(o.begin),
                                       ExactPosition(o.end))
            record.features.append(
                SeqFeature(location,
                    type='mRNA',
                    strand=o.strand,
                    qualifiers=dict(
                        note='putative, confidence %d%%' % o.confidence,
                        operon='rnas-%d' % i)))

        record.features.sort(key=lambda f: f.location.start.position)

        xgb_file = open(self.analysis.xgenbankfile_path, "w")

        SeqIO.write(record, xgb_file, "genbank")
    #
    ############################################################################

    __iter__ = lambda self: self._stages

WorkerStages._stages = _stages
del _stages

class UnknownInputfileTypeException(Exception):
    pass
