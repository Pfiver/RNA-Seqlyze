from os.path import join, exists, isdir

from .utils import OrderedClass

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

    __metaclass__ = OrderedClass

    def determine_inputfile_type(self):
        def getit(header):
            if header[0] == '@': return 'fastq'
            elif header[:8] == 'NCBI.sra': return 'sra'
            else: raise UnknownInputfileTypeException()
        self.analysis.inputfile_type = \
                getit(open(self.analysis.inputfile_path).read(8))
        self.session.commit()

    def convert_input_file(self):
        fq_path = self.analysis.inputfile_fqpath
        if not exists(fq_path):
            import os
            os.chdir(self.input_data_dir)
            if self.analysis.inputfile_name:
                sra_name = self.analysis.inputfile_name
            else:
                sra_name = self.analysis.rnaseq_run.sra_name
            cmd = "fastq-dump", sra_name
            self.log.info("converting %s" % sra_name)
            from subprocess import Popen, PIPE
            proc = Popen(cmd)#, stdout=PIPE, stderr=PIPE)
            out, err = proc.communicate()
            if proc.returncode != 0:
                raise Exception("%s failed" % (cmd,))

    def fetch_gb(self):
        from os import path
        acc = self.analysis.org_accession
        self.analysis.create_gb_data_dir()
        out_path = join(self.gb_data_dir, acc + ".gb")
        if not exists(out_path):
            self.log.info("Fetching '%s' from entrez..." % acc)
            # TODO: do this earlier
            gb_id = efetch.get_nc_id(acc)
            efetch.fetch_nc_gb(gb_id, open(out_path, "w"))
            self.log.info("...done")

    def gb2fasta(self):
        from os import path
        acc = self.analysis.org_accession
        gb_path = join(self.gb_data_dir, acc + ".gb")
        fa_path = join(self.gb_data_dir, acc + ".fa")
        if not exists(fa_path):
            self.log.info("Converting '%s' to fasta format..." % acc)
            import Bio.SeqIO
            record = Bio.SeqIO.parse(open(gb_path), "genbank").next()
            record.id = "chr" # required (!) for ucsc browser
            Bio.SeqIO.write(record, open(fa_path, "w"), "fasta")

    def bowtie_build(self):
        from os import path
        acc = self.analysis.org_accession
        bt2_path = join(self.gb_data_dir, acc + ".1.bt2")
        if not exists(bt2_path):
            import os
            os.chdir(self.gb_data_dir)
            cmd = "bowtie2-build", acc + ".fa", acc
            self.log.info("bowtie2-build %s" % acc)
            from subprocess import Popen, PIPE
            proc = Popen(cmd)#, stdout=PIPE, stderr=PIPE)
            out, err = proc.communicate()
            if proc.returncode != 0:
                raise Exception("%s failed" % (cmd,))

    def tophat(self):
        from os import path
        acc = self.analysis.org_accession
        if not isdir(join(self.data_dir, "tophat-output")):
            import os
            os.chdir(self.data_dir)
            n_cpus = os.sysconf("SC_NPROCESSORS_ONLN")
            acc_path = join(self.gb_data_dir, acc)
            fq_name = self.analysis.inputfile_fqname
            cmd = "tophat", "-p", str(n_cpus), \
                    "-o", "tophat-output", acc_path, fq_name
            from subprocess import Popen, PIPE
            proc = Popen(cmd)#, stdout=PIPE, stderr=PIPE)
            out, err = proc.communicate()
            if proc.returncode != 0:
                raise Exception("%s failed" % (cmd,))

    def galaxy_upload(self):
        if self.analysis.galaxy_bam_id:
            return
        from os import path
        from rnaseqlyze import galaxy
        bam_path = join(self.data_dir, 
                             "tophat-output", "accepted_hits.bam")
        if self.analysis.inputfile_name: # uploaded input file
            srr_name = self.analysis.inputfile_name.rsplit(".", 1)[0]
        else: # input is an SRR Indentifier
            srr_name = self.analysis.rnaseq_run.srr
        # FIXME: names are not unique on galaxy:
        # is "%s_%s" % (srr_name, self.analysis.org_accession) good enough ?
        galaxy_bam_name = "%s_%s" % (srr_name, self.analysis.org_accession)
        bam_file = open(bam_path)
        self.log.info("uploading %s to galaxy server %s ..." % (
                                 bam_path,           galaxy.hostname))
        self.analysis.galaxy_bam_id = galaxy.upload(bam_file, galaxy_bam_name)
        self.log.info("done.")
        bam_file.close()
        self.session.commit()

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

        operons = dict(beg=0, end=100, strand=1),

        from os import path
        acc = self.analysis.org_accession
        gb_path = join(self.gb_data_dir, acc + ".gb")

        from Bio import SeqIO
        from Bio.SeqFeature import \
                SeqFeature, FeatureLocation

        record = SeqIO.parse(open(gb_path), "genbank").next()

        for oper in operons:
            location = FeatureLocation(ExactPosition(oper.beg),
                                       ExactPosition(oper.end))
            record.features.append(
                SeqFeature(location, type='mRNA', strand=oper.strand))

        record.features.sort(key=lambda f: f.location.start.position)
        SeqIO.write(record, open(xgb_path, "w"), "genbank")


class UnknownInputfileTypeException(Exception):
    pass
