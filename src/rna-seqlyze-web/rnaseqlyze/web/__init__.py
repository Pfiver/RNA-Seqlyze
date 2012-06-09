from sqlalchemy import create_engine
from pyramid.config import Configurator

import rnaseqlyze
from rnaseqlyze.core.orm import DBSession


def main(global_config, **settings):
    """
    Return a Pyramid(!) WSGI application.
    """
    engine = create_engine(rnaseqlyze.db_url)
    DBSession.configure(bind=engine)

    config = Configurator(settings=settings)

    config.add_subscriber('rnaseqlyze.web.add_base_template',
                          'pyramid.events.BeforeRender')

    config.add_route('home', '/')
    config.add_route('analysis', '/analysis')
    for path in 'img', 'css', 'js':
        config.add_static_view(path, path)
    config.scan()

    return config.make_wsgi_app()



from pyramid.renderers import get_renderer

def add_base_template(event):
    base = get_renderer('templates/base.pt').implementation()
    event.update({'base': base})
