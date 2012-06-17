import logging
log = logging.getLogger(__name__)

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker, scoped_session

DBSession = scoped_session(sessionmaker())

import rnaseqlyze
from ..core import service
from ..core.orm import Analysis


# Allways commit the DBSession after making changes!
# Otherwise the db will stay locked by this thread and
# other threads/processes trying to access it will block

class Manager(object):
    def __init__(self):
        engine = create_engine(rnaseqlyze.db_url)
        DBSession.configure(bind=engine)
        self.worker = None

    def get_analysis(self, id):
        # refresh session
        DBSession.commit()
        return DBSession.query(Analysis).get(id)

    def start_analysis(self, analysis):
        if analysis.started:
            raise AnalysisAlreadyStartedException(analysis.id)
        analysis.started = True
        DBSession.commit()
        if not self.worker:
            self.worker = AnalysisOrchestrator(analysis.id)
            self.worker.start()

class AnalysisAlreadyStartedException(Exception):
    pass


from threading import Thread
class AnalysisOrchestrator(Thread):

    def __init__(self, analysis_id):
        Thread.__init__(self)
        self.analysis_id = analysis_id

    def determine_inputfile_type(self):
        def getit(header):
            if header[0] == '@': return 'fastq'
            elif header[:8] == 'NCBI.sra': return 'sra'
            else: raise UnknownInputfileTypeException()
        path = service.get_inputfile_path(self.analysis)
        self.analysis.inputfile_type = getit(open(path).read(8))
        DBSession.commit()

    def run(self):
        self.analysis = DBSession.query(Analysis).get(self.analysis_id)
        self.determine_inputfile_type()

class UnknownInputfileTypeException(Exception):
    pass
