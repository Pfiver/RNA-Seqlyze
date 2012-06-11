from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker, scoped_session

import rnaseqlyze

DBSession = scoped_session(sessionmaker())

def main():
    engine = create_engine(rnaseqlyze.db_url)
    DBSession.configure(bind=engine)
