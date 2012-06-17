import logging
log = logging.getLogger(__name__)
from threading import Thread

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

import rnaseqlyze
from ..core import service
from ..core.orm import Analysis


# Allways commit the DBSession after making changes!
# Otherwise the db will stay locked by this thread and
# other threads/processes trying to access it will block

DBSession = sessionmaker(create_engine(rnaseqlyze.db_url))

class Manager(object):
    def __init__(self):
        self.worker = Thread()
        self.session = DBSession()

    def get_analysis(self, id):
        return self.session.query(Analysis).get(id)

    def start_analysis(self, analysis):
        if analysis.started:
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
        self._thread_init()
        self._determine_inputfile_type()

    def _thread_init(self):
        self.session = DBSession()
        self.analysis = self.session.query(Analysis).get(self.analysis_id)
        self.analysis.started = True
        self.session.commit()

    def _determine_inputfile_type(self):
        def getit(header):
            if header[0] == '@': return 'fastq'
            elif header[:8] == 'NCBI.sra': return 'sra'
            else: raise UnknownInputfileTypeException()
        path = service.get_inputfile_path(self.analysis)
        self.analysis.inputfile_type = getit(open(path).read(8))
        self.session.commit()

class UnknownInputfileTypeException(Exception):
    pass
