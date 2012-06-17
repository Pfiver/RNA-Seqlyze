import logging
log = logging.getLogger(__name__)

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

DBSession = sessionmaker()

import rnaseqlyze
from ..core.orm import Analysis


class Manager(object):

    def __init__(self):
        engine = create_engine(rnaseqlyze.db_url)
        DBSession.configure(bind=engine)
        self.db_session = DBSession()

    def get_analysis(self, id):
        # refresh session
        self.db_session.commit()
        return self.db_session.query(Analysis).get(id)

    def start_analysis(self, analysis):
        if analysis.started:
            raise AnalysisAlreadyStartedException(analysis.id)
        analysis.started = True
        self.db_session.commit()


class AnalysisAlreadyStartedException(Exception):
    pass
