from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

import rnaseqlyze

DBSession = sessionmaker(create_engine(rnaseqlyze.db_url))

def main():
    session = DBSession()
