"""
RNA-Seqlyze Worker Daemon Core

Worker parent class with basic infrastructure to run the
various analysis steps defined in :class:`~.WorkerStages`.
"""

import logging
log = logging.getLogger(__name__)

from threading import Thread
from logging import Formatter
from logging import StreamHandler
from StringIO import StringIO
from contextlib import contextmanager

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

import rnaseqlyze
from rnaseqlyze import efetch
from rnaseqlyze.core.entities import Analysis, StageLog
from rnaseqlyze.worker.stages import WorkerStages

DBSession = sessionmaker()

log_format = "%(levelname)-5.5s [%(name)s] %(message)s"

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

class StageLogStream(object):
    def __init__(self, analysis, stage, session):
        self.stage_log = StageLog(analysis=analysis, stage=stage, text="")
        self.session = session
        session.add(self.stage_log)
        session.commit()
    def write(self, data):
        self.stage_log.text += data
        self.session.commit()

class AnalysisAlreadyStartedException(Exception):
    pass

class ManagerBusyException(Exception):
    pass

class Worker(Thread, WorkerStages):
    """
    The Worker
    """

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
        self.log_handler = StreamHandler(self.logfile)
        self.log_handler.setFormatter(Formatter(log_format))
        self.log.addHandler(self.log_handler)
        self.log.info("starting work on analysis #%d" % self.analysis_id)
        self.analysis.finished = False
        self.analysis.started = True
        self.session.commit()

    @contextmanager
    def _stage_log_manager(self, stage):
        handler = StreamHandler(StageLogStream(
                                        self.analysis, stage, self.session))
        self.log.addHandler(handler)
        yield
        self.log.removeHandler(handler)

    def run(self):

        # TODO: invent a way to avoid calling stages that won't do anything
        #       maybe @stages -> @stages(condition) something ...

        self._thread_init()
        self.analysis.error = None
        self.session.commit()
        try:
            for stage in WorkerStages:
                self.log.info("=== %s ===" % stage.func_name)
                with self._stage_log_manager(stage.func_name):
                    self.analysis.stage = stage.func_name
                    self.session.commit()
                    stage(self)
        except Exception, e:
            self.analysis.error = repr(e)
            raise
        finally:
            if self.analysis.error:
                self.log.error("error processing analysis #%d: %r" % (
                                                           self.analysis.id,
                                                           self.analysis.error))
            else:
                self.log.info("work on analysis #%d completed" % (
                                                 self.analysis.id,))
            self.analysis.finished = True
            self.session.commit()
            self.log.removeHandler(self.log_handler)
            self.logfile.close()
