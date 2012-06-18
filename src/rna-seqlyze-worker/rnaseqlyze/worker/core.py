import logging
log = logging.getLogger(__name__)
from threading import Thread

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

import rnaseqlyze
from rnaseqlyze import efetch
from rnaseqlyze.core import service
from rnaseqlyze.core.orm import Analysis


# Allways commit the DBSession after making changes!
# Otherwise the db will stay locked by this thread and
# other threads/processes trying to access it will block

DBSession = sessionmaker(create_engine(rnaseqlyze.db_url))

class Manager(object):
    def __init__(self):
        self.worker = Thread()
        self.session = DBSession()

    def get_analysis(self, id):
        obj = self.session.query(Analysis).get(id)
        self.session.refresh(obj)
        return obj

    def start_analysis(self, analysis, re=False):
        if analysis.started and not re:
            raise AnalysisAlreadyStartedException(analysis)
        if self.worker.is_alive():
            raise ManagerBusyException

        self.worker = AnalysisOrchestrator(analysis.id)
        self.worker.start()

class AnalysisAlreadyStartedException(Exception):
    pass

class ManagerBusyException(Exception):
    pass


class AnalysisOrchestrator(Thread):

    def __init__(self, analysis_id):
        Thread.__init__(self)
        self.analysis_id = analysis_id

    def run(self):
        log.info("starting work on analysis #%d" % self.analysis_id)
        self._thread_init()
        self._determine_inputfile_type()
        self._get_input_header()
        self._fetch_gb()
        self._gb2fasta()
        self._bowtie_build()

    def _thread_init(self):
        self.session = DBSession()
        self.analysis = self.session.query(Analysis).get(self.analysis_id)
        self.analysis.started = True
        self.session.commit()
        self.data_dir = service.get_data_dir(self.analysis)
        self.shared_data_dir = service.get_shared_data_dir(self.analysis)

    def _determine_inputfile_type(self):
        def getit(header):
            if header[0] == '@': return 'fastq'
            elif header[:8] == 'NCBI.sra': return 'sra'
            else: raise UnknownInputfileTypeException()
        path = service.get_inputfile_path(self.analysis)
        self.analysis.inputfile_type = getit(open(path).read(8))
        self.session.commit()

    def _get_input_header(self):
        type = self.analysis.inputfile_type
        path = service.get_inputfile_path(self.analysis)
        if type == 'fastq':
            f = open(path)
            data = [f.readline() for i in range(4)]
        elif type == 'sra':
            cmd = "fastq-dump", "--stdout", path
            log.info("converting %s" % path)
            from subprocess import Popen, PIPE
            proc = Popen(cmd, stdout=PIPE)
            data = [proc.stdout.readline() for i in range(4)]
            proc.kill()
            proc.wait()
        else:
            raise
        log.debug("got header: %s ..." % data[0][:-1])
        self.analysis.inputfile_header = "".join(data)
        self.session.commit()

    def _fetch_gb(self):
        from os import path
        acc = self.analysis.org_accession
        out_path = path.join(self.shared_data_dir, acc + ".gb")
        if not path.exists(out_path):
            log.info("Fetching '%s' from entrez..." % acc)
            # TODO: do this earlier
            gb_id = efetch.get_nc_id(acc)
            efetch.fetch_nc_gb(gb_id, open(out_path, "w"))
            log.info("...done")

    def _gb2fasta(self):
        from os import path
        acc = self.analysis.org_accession
        gb_path = path.join(self.shared_data_dir, acc + ".gb")
        fa_path = path.join(self.shared_data_dir, acc + ".fa")
        if not path.exists(fa_path):
            log.info("Converting '%s' to fasta format..." % acc)
            import Bio.SeqIO
            parse = Bio.SeqIO.parse(open(gb_path), "genbank")
            Bio.SeqIO.write(parse, open(fa_path, "w"), "fasta") 

    def _bowtie_build(self):
        from os import path
        acc = self.analysis.org_accession
        bt2_path = path.join(self.shared_data_dir, acc + ".1.bt2")
        if not path.exists(bt2_path):
            import os
            os.chdir(self.shared_data_dir)
            cmd = "bowtie2-build", acc + ".fa", acc
            log.info("bowtie2-build %s" % acc)
            from subprocess import Popen, PIPE
            proc = Popen(cmd, stdout=PIPE)
            out, err = proc.communicate()
            if proc.returncode != 0:
                raise Exception("%s failed" % (cmd,))

class UnknownInputfileTypeException(Exception):
    pass
