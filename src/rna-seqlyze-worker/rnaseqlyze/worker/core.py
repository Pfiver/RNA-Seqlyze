"""
RNA-Seqlyze Worker Daemon Core

    -- this is where things are actually getting done! :-)
"""

import logging
log = logging.getLogger(__name__)

from threading import Thread
from logging import Formatter
from logging import StreamHandler

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

import rnaseqlyze
from rnaseqlyze import efetch
from .stages import WorkerStages
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

class Worker(Thread, WorkerStages):
    """
    The Worker
    """

    # Allways commit the DBSession after making changes!
    # Otherwise the db will stay locked by this thread and
    # other threads/processes trying to access it will block

    def __init__(self, analysis):
        Thread.__init__(self)
        self.analysis_id = analysis.id

    def _thread_init(self):
        from os import path
        self.session = DBSession(bind=create_engine(rnaseqlyze.db_url))
        self.analysis = self.session.query(Analysis).get(self.analysis_id)

        self.logfile = open(path.join(
                self.analysis.data_dir, "rna-seqlyze-worker.log"), "w")
        self.log = logging.getLogger("%s.Worker(Analysis #%d)" % (
                                      __name__,           self.analysis_id))
        h = StreamHandler(self.logfile)
        h.setFormatter(Formatter("%(levelname)-5.5s [%(name)s] %(message)s"))
        self.log.addHandler(h)
        self.log.info("starting work on analysis #%d" % self.analysis_id)
        self.analysis.started = True
        self.session.commit()

        self.data_dir = self.analysis.data_dir
        self.gb_data_dir = self.analysis.gb_data_dir
        self.input_data_dir = self.analysis.input_data_dir

    def run(self):
        self._thread_init()
        try:
            for stage in WorkerStages.members:
                self.log.info("=== %s ===" % stage.func_name)
                stage(self)
        except Exception, e:
            self.analysis.error = str(e)
            raise
        finally:
            self.log.info("work on analysis #%d completed" % self.analysis_id)
            self.analysis.finished = True
            self.session.commit()
