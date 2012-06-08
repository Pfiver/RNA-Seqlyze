from sqlalchemy import create_engine
from pyramid.config import Configurator

import rnaseqlyze
from rnaseqlyze.orm import DBSession

def main(global_config, **settings):
    """
    Return a Pyramid(!) WSGI application.
    """

    engine = create_engine(rnaseqlyze.db_url)
    DBSession.configure(bind=engine)

    config = Configurator(settings=settings)
    config.scan()

    return config.make_wsgi_app()
