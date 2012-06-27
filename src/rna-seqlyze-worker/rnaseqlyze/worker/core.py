import logging
log = logging.getLogger(__name__)

from threading import Thread

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

import rnaseqlyze
from rnaseqlyze import efetch
from rnaseqlyze.core.orm import Analysis

DBSession = sessionmaker()

class Manager(object):
    def __init__(self):
        self.worker = Thread()

    def analysis_requested(self, analysis, re=False):
        if analysis.started and not re:
            raise AnalysisAlreadyStartedException
        if self.worker.is_alive():
            raise ManagerBusyException
        self.worker = Worker(analysis)
        self.worker.start()

class AnalysisAlreadyStartedException(Exception):
    pass

class ManagerBusyException(Exception):
    pass


class Worker(Thread):

    # Allways commit the DBSession after making changes!
    # Otherwise the db will stay locked by this thread and
    # other threads/processes trying to access it will block

    def __init__(self, analysis):
        Thread.__init__(self)
        self.analysis_id = analysis.id

    def _thread_init(self):
        log.info("starting work on analysis #%d" % self.analysis_id)
        self.session = DBSession(bind=create_engine(rnaseqlyze.db_url))
        self.analysis = self.session.query(Analysis).get(self.analysis_id)
        self.analysis.started = True
        self.session.commit()
        self.data_dir = self.analysis.data_dir
        self.gb_data_dir = self.analysis.gb_data_dir
        self.input_data_dir = self.analysis.input_data_dir

    def run(self):
        self._thread_init()
        self._determine_inputfile_type()
        self._convert_input_file()
        self._fetch_gb()
        self._gb2fasta()
        self._bowtie_build()
        self._tophat()
        self._galaxy_upload()

    def _determine_inputfile_type(self):
        def getit(header):
            if header[0] == '@': return 'fastq'
            elif header[:8] == 'NCBI.sra': return 'sra'
            else: raise UnknownInputfileTypeException()
        self.analysis.inputfile_type = \
                getit(open(self.analysis.inputfile_path).read(8))
        self.session.commit()

    def _convert_input_file(self):
        from os import path
        fq_path = self.analysis.inputfile_fqpath
        if not path.exists(fq_path):
            import os
            os.chdir(self.input_data_dir)
            if self.analysis.inputfile_name:
                sra_name = self.analysis.inputfile_name
            else:
                sra_name = self.analysis.rnaseq_run.srr + '.sra'
            cmd = "fastq-dump", sra_name
            log.info("converting %s" % sra_name)
            from subprocess import Popen, PIPE
            proc = Popen(cmd)#, stdout=PIPE, stderr=PIPE)
            out, err = proc.communicate()
            if proc.returncode != 0:
                raise Exception("%s failed" % (cmd,))

    def _fetch_gb(self):
        from os import path
        acc = self.analysis.org_accession
        self.analysis.create_gb_data_dir()
        out_path = path.join(self.gb_data_dir, acc + ".gb")
        if not path.exists(out_path):
            log.info("Fetching '%s' from entrez..." % acc)
            # TODO: do this earlier
            gb_id = efetch.get_nc_id(acc)
            efetch.fetch_nc_gb(gb_id, open(out_path, "w"))
            log.info("...done")

    def _gb2fasta(self):
        from os import path
        acc = self.analysis.org_accession
        gb_path = path.join(self.gb_data_dir, acc + ".gb")
        fa_path = path.join(self.gb_data_dir, acc + ".fa")
        if not path.exists(fa_path):
            log.info("Converting '%s' to fasta format..." % acc)
            import Bio.SeqIO
            record = Bio.SeqIO.parse(open(gb_path), "genbank").next()
            record.id = "chr" # required (!) for ucsc browser
            Bio.SeqIO.write(record, open(fa_path, "w"), "fasta")

    def _bowtie_build(self):
        from os import path
        acc = self.analysis.org_accession
        bt2_path = path.join(self.gb_data_dir, acc + ".1.bt2")
        if not path.exists(bt2_path):
            import os
            os.chdir(self.gb_data_dir)
            cmd = "bowtie2-build", acc + ".fa", acc
            log.info("bowtie2-build %s" % acc)
            from subprocess import Popen, PIPE
            proc = Popen(cmd)#, stdout=PIPE, stderr=PIPE)
            out, err = proc.communicate()
            if proc.returncode != 0:
                raise Exception("%s failed" % (cmd,))

    def _tophat(self):
        from os import path
        acc = self.analysis.org_accession
        if not path.isdir(path.join(self.data_dir, "tophat-output")):
            import os
            os.chdir(self.data_dir)
            n_cpus = os.sysconf("SC_NPROCESSORS_ONLN")
            acc_path = path.join(self.gb_data_dir, acc)
            fq_name = self.analysis.inputfile_fqname
            cmd = "tophat", "-p", str(n_cpus), \
                    "-o", "tophat-output", acc_path, fq_name
            from subprocess import Popen, PIPE
            proc = Popen(cmd)#, stdout=PIPE, stderr=PIPE)
            out, err = proc.communicate()
            if proc.returncode != 0:
                raise Exception("%s failed" % (cmd,))

    def _galaxy_upload(self):
        if self.analysis.galaxy_bam_id:
            return
        from os import path
        from rnaseqlyze import galaxy
        bam_path = path.join(self.data_dir, 
                             "tophat-output", "accepted_hits.bam")
        srr_name = self.analysis.inputfile_name.rsplit(".", 1)[0]
        # FIXME: names are not unique on galaxy:
        # is "%s_%s" % (srr_name, self.analysis.org_accession) good enough ?
        galaxy_bam_name = "%s_%s" % (srr_name, self.analysis.org_accession)
        bam_file = open(bam_path)
        log.info("uploading %s to galaxy server %s ..." % (
                            bam_path,           galaxy.hostname))
        self.analysis.galaxy_bam_id = galaxy.upload(bam_file, galaxy_bam_name)
        log.info("done.")
        bam_file.close()
        self.session.commit()


class UnknownInputfileTypeException(Exception):
    pass
