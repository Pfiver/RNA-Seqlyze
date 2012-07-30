"""
RNA-Seqlyze Worker Daemon Core

Worker parent class with basic infrastructure to run the
various analysis steps defined in :class:`~.WorkerStages`.
"""

import logging
log = logging.getLogger(__name__)
root_logger = logging.getLogger()

from threading import Thread
from logging import Formatter
from logging import StreamHandler
from StringIO import StringIO
from contextlib import contextmanager

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

import rnaseqlyze
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
        self.log_handler = StreamHandler(self.logfile)
        self.log_handler.setFormatter(Formatter(log_format))
        root_logger.addHandler(self.log_handler)
        log.info("starting work on analysis #%d" % self.analysis_id)
        self.analysis.finished = False
        self.analysis.started = True
        self.analysis.error = None
        self.session.commit()

    @contextmanager
    def _stage_log_manager(self, stage):
        handler = StreamHandler(
                    StageLogStream(self.analysis, stage, self.session))
        root_logger.addHandler(handler)
        try:
            yield
        finally:
            root_logger.removeHandler(handler)

    def run(self):

        # TODO: invent a way to avoid calling stages that won't do anything
        #       maybe @stages -> @stages(condition) something ...

        self._thread_init()
        try:
            for stage in self.stages:
                if not stage.should_run(self):
                    continue
                log.info("=== %s ===" % stage.func_name)
                with self._stage_log_manager(stage.func_name):
                    self.analysis.stage = stage.func_name
                    self.session.commit()
                    stage(self)
        except Exception, e:
            self.analysis.error = repr(e)
            raise
        finally:
            if self.analysis.error:
                log.error(self.analysis.error)
            else:
                log.info("analysis finished")

            self.analysis.finished = True
            self.session.commit()
            root_logger.removeHandler(self.log_handler)
            self.logfile.close()
