from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker, scoped_session
from pyramid.config import Configurator

import rnaseqlyze

DBSession = scoped_session(sessionmaker())

def main():
    engine = create_engine(rnaseqlyze.db_url)
    DBSession.configure(bind=engine)

def main(global_config, **settings):
    """
    Return a Pyramid(!) WSGI application.
    """
    engine = create_engine(rnaseqlyze.db_url)
    DBSession.configure(bind=engine)

    config = Configurator(settings=settings)
    config.scan()

    return config.make_wsgi_app()
