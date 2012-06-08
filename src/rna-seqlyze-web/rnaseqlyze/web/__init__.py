from pyramid.config import Configurator
from sqlalchemy import engine_from_config

from rnaseqlyze.core.orm import DBSession

def main(global_config, **settings):
    """
    Returns a Pyramid WSGI application.
    """
    engine = engine_from_config(settings, 'sqlalchemy.')
    DBSession.configure(bind=engine)

    config = Configurator(settings=settings)
    config.add_route('base', '/')
    for path in 'img', 'css', 'js':
        config.add_static_view(path, path)
    config.scan()
    return config.make_wsgi_app()

