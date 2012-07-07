"""
RNA-Seqlyze Worker Stages

    -- **this** is where things are actually getting done! :-)
"""
import logging
log = logging.getLogger(__name__)

import os
from os.path import join, exists, isdir, relpath
from subprocess import Popen, PIPE
from StringIO import StringIO
from threading import Thread
from urllib import quote

import pysam

from Bio import SeqIO
from Bio.SeqFeature import \
        SeqFeature, FeatureLocation, ExactPosition

from psutil import cpu_percent

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
_stage_conds = {}
def stage(method):
    """
    Just a small helper to collect the stages in the order defined.

    To add a new stage, simply add a method to :class:`~WorkerStages`.
    It will be automatically executed for all new analyses.
    """
    if method.func_name in _stage_conds:
        method.should_run = _stage_conds[method.func_name]
        del _stage_conds[method.func_name]
    else:
        method.should_run = lambda self: True
    _stages.append(method)
    return method

def stage_cond(method):
    """
    Stage Condition
    
     - must be declared before @stage
     - must return true for the @stage of the same name to run
    """
    _stage_conds[method.func_name] = method

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
        # can't wait() on subprocess with a timeout, alas start up
        # a 2nd thread to do it and join() on that one with a timeout
        log.info("forking subprocess: $ %s" % ' '.join(map(repr, cmd)))
        proc = Popen(cmd, stdout=self.logfile, stderr=self.logfile)
        waiter = Thread(target=proc.wait)
        cpu_percent(0, True)
        waiter.start()
        waiter.join(15)
        while waiter.is_alive():
            log.info("subprocess still running - system load: " +
                     " / ".join(("%d%%" % p for p in cpu_percent(0, True))))
            waiter.join(15)
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

    @stage_cond
    def determine_inputfile_type(self):
        return self.analysis.inputfile_uploaded
    @stage
    def determine_inputfile_type(self):
        _8bytes = open(self.analysis.inputfile_path).read(8)
        log.info("first 8 bytes of input data: %r" % _8bytes)
        self.analysis.inputfile_type = (
                     'fastq' if _8bytes[0] == '@'
                else 'sra'   if header == 'NCBI.sra' else None)
        self.session.commit()
        if not self.analysis.inputfile_type:
            raise Exception("Unknown input data type")

    @stage_cond
    def fetch_srr(self):
        # don't download if private
        # file uploaded or srr already in cache
        return not self.analysis.inputfile_uploaded \
               and not os.path.exists(self.analysis.rnaseq_run.sra_path)
    @stage
    def fetch_srr(self):
        self.analysis.rnaseq_run.download()

    @stage_cond
    def convert_input_file(self):
        return not exists(self.analysis.inputfile_fq_path)
    @stage
    def convert_input_file(self):
        os.chdir(self.analysis.input_data_dir)
        self.log_cmd("fastq-dump", self.analysis.inputfile_name)
        log.debug("created %s" % self.analysis.inputfile_fq_path)

    @stage_cond
    def fetch_genbank_file(self):
        return not exists(self.analysis.genbankfile_path)
    @stage
    def fetch_genbank_file(self):
        if not exists(self.analysis.genbank_data_dir):
            os.makedirs(self.analysis.genbank_data_dir)
        log.info("Fetching '%s' from entrez..." %
                                 self.analysis.org_accession)
        gb_id = efetch.get_nc_id(self.analysis.org_accession)
        efetch.fetch_nc_gb(gb_id, open(self.analysis.genbankfile_path, "w"))
        log.info("...done")

    @stage
    def read_genbank_file(self):
        self.genbank_record = SeqIO.parse(open(self.analysis \
                                    .genbankfile_path), "genbank").next()
        ngenes = sum(1 for f in self.genbank_record.features
                       if f.type == 'gene')
        log.info("genbank file lists %d genes" % ngenes)

    @stage_cond
    def genbank_to_fasta(self):
        return not exists(self.analysis.genbankfile_fa_path)
    @stage
    def genbank_to_fasta(self):
        log.info("Converting '%s' to fasta format" %
                                           self.analysis.genbankfile_name)
        record = self.genbank_record
        saved_id = record.id
        record.id = "chr" # make ucsc browser custom tracks work
        SeqIO.write(record, open(
            self.analysis.genbankfile_fa_path, "w"), "fasta")
        record.id = saved_id

    @stage_cond
    def bowtie_build(self):
        return not exists(join(self.analysis.genbank_data_dir,
                               self.analysis.genbankfile_base_name + ".1.bt2"))
    @stage
    def bowtie_build(self):
        os.chdir(self.analysis.genbank_data_dir)
        self.log_cmd("bowtie2-build", self.analysis.genbankfile_fa_name,
                                      self.analysis.genbankfile_base_name)

    @stage_cond
    def tophat(self):
        return not exists(join(self.analysis.data_dir,
                               "tophat-output", "accepted_hits.bam"))
    @stage
    def tophat(self):
        os.chdir(self.analysis.data_dir)
        n_cpus = os.sysconf("SC_NPROCESSORS_ONLN")
        fq = relpath(self.analysis.inputfile_fq_path)
        gb = relpath(join(self.analysis.genbank_data_dir,
                          self.analysis.genbankfile_base_name))
        self.log_cmd("tophat", "-p", str(n_cpus), "-o", "tophat-output", gb, fq)

    @stage_cond
    def create_coverage_track(self):
        return not exists(join(self.analysis.data_dir, "coverage.bigwig"))
    @stage
    def create_coverage_track(self):
        os.chdir(self.analysis.data_dir)
        # the script automatically converts it's
        # output to bigwig if it finds kent's wigToBigWig
        self.log_cmd("bam_to_wiggle.py", "-o", "coverage.bigwig",
                                    "tophat-output/accepted_hits.bam")

    @stage_cond
    def genbank_to_ptt(self):
        return not exists(join(self.analysis.genbank_data_dir,
                               self.genbank_record.id + ".ptt"))
    @stage
    def genbank_to_ptt(self):
        ptt_name = self.genbank_record.id + ".ptt"
        ptt_path = join(self.analysis.genbank_data_dir, ptt_name)
        os.symlink(ptt_name, join(self.analysis.genbank_data_dir, "chr.ptt"))
        log.debug("converting %s to ptt" % self.analysis.genbankfile_name)
        ptt_file = open(ptt_path, "w")
        gb_file = open(self.analysis.genbankfile_path)
        gb2ptt.gb2ptt(gb_file, ptt_file)
        ptt_file.close()
        gb_file.close()

    @stage_cond
    def transterm_hp(self):
        return not exists(join(self.analysis.data_dir,
                               "hp_terminators.bigbed"))
    @stage
    def transterm_hp(self):
        os.chdir(self.analysis.data_dir)
        log.debug("running transterm")

        tt_out = open("transterm_hp.out", "w+")
        # --min-conf=n n is the cut-off confidence value,
        #              between 0 and 100, the default is 76
        tt_args = ("--min-conf=47",
                   self.analysis.genbankfile_fa_path,
                   relpath(join(self.analysis.genbank_data_dir, "chr.ptt")))
        transterm.run(tt_args, out=tt_out, err=self.logfile)
        tt_out.seek(0)
        # keep a copy in memory
        self.hp_terminators = list(transterm.iterator(tt_out))
        tt_out.seek(0)
        log.info("found {0} possible hairpin terminators"
                                 .format(len(self.hp_terminators)))
        # create a bed track
        bed_file = open("hp_terminators.bed", "w")
        transterm.tt2bed(tt_out, bed_file)
        bed_file.close()
        tt_out.close()

        log.debug("running bedToBigBed")

        # convert it to bigBed
        chrs = open("chrom.sizes", "w")
        chrs.write("chr %d" % len(self.genbank_record.seq))
        chrs.close()
        self.log_cmd("bedToBigBed", "hp_terminators.bed",
                         "chrom.sizes", "hp_terminators.bigbed")

    @stage
    def predict_operons(self):

        # extract the coverage data from the bam track created by tophat
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

        log.debug("maximum coverage: %d" % self.max)
        log.debug("number of bases covered by short reads: %d/%d" % (
                                    self.covered, len(self.genbank_record.seq)))

        # available objects at this point
        # -------------------------------
        #
        # - self.genbank_record: Biopython SeqIO.parse()d genbank file
        #
        # - self.coverage: [n,n,n,n,...] / len = len(self.genbank_record.seq)
        #
        # - self.max: max(n)
        #
        # - self.hpterminators: ((id, begin, end, strand, confidence), ...)
        #                         str, str, str, str (1/-), int
        #

        # FIXME: do some magic here
        self.operons = Operon(begin=0, end=100, strand=1, confidence=10),
        self.operons = Operon(begin=200, end=300, strand=1, confidence=50),
        self.operons = Operon(begin=400, end=500, strand=1, confidence=100),

        # create a bed track
        track_name = "rna-seqlyze-operon_predictions"
        os.chdir(self.analysis.data_dir)
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

        # convert it to bigBed
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
            log.info("uploading accepted_hits.bam to galaxy")
            self.analysis.galaxy_bam = GalaxyDataset(
                id=galaxy.upload(open(bam_path), self.bam_name))
            log.info("...done - id: %s" % self.analysis.galaxy_bam.id)
            self.session.commit()

        if not self.analysis.galaxy_coverage:
            coverage_path = join(self.analysis.data_dir, "coverage.bigwig")
            log.info("uploading coverage.bigwig to galaxy")
            self.analysis.galaxy_coverage = GalaxyDataset(
                id=galaxy.upload(open(coverage_path), self.coverage_name))
            log.info("...done - id: %s" % self.analysis.galaxy_coverage.id)
            self.session.commit()

        if not self.analysis.galaxy_hp_terms:
            hp_terms_path = join(self.analysis.data_dir,
                                 "hp_terminators.bigbed")
            log.info("uploading hp_terminators.bigbed to galaxy")
            self.analysis.galaxy_hp_terms = GalaxyDataset(
                id=galaxy.upload(open(hp_terms_path), self.hp_terms_name))
            log.info("...done - id: %s" % self.analysis.galaxy_hp_terms.id)
            self.session.commit()

        if not self.analysis.galaxy_pr_operons:
            track_filename = "rna-seqlyze-operon_predictions.bigbed"
            pr_operons_path = join(self.analysis.data_dir, track_filename)
            log.info("uploading %s to galaxy" % track_filename)
            self.analysis.galaxy_pr_operons = GalaxyDataset(
                id=galaxy.upload(open(pr_operons_path), self.pr_operons_name))
            log.info("...done - id: %s" %
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

        # FIXME: this cries for refactoring -- with logging!

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

        log.info("augmenting genbank file %s with putative operons" %
                                          self.analysis.genbankfile_name)

        for i, o in enumerate(self.operons):
            location = FeatureLocation(ExactPosition(o.begin),
                                       ExactPosition(o.end))
            self.genbank_record.features.append(
                SeqFeature(location,
                    type='mRNA',
                    strand=o.strand,
                    qualifiers=dict(
                        note='putative, confidence %d%%' % o.confidence,
                        operon='rnas-%d' % i)))

        self.genbank_record.features.sort(
                key=lambda f: f.location.start.position)

        xgb_file = open(self.analysis.xgenbankfile_path, "w")

        SeqIO.write(self.genbank_record, xgb_file, "genbank")
    #
    ############################################################################

assert not _stage_conds, "@stage_cond's must be declared before @stage's"
WorkerStages.stages = _stages
del _stages, _stage_conds
