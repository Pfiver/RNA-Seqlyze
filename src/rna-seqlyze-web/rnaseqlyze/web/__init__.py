from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker, scoped_session
from pyramid.config import Configurator
from zope.sqlalchemy import ZopeTransactionExtension

import rnaseqlyze

# why ZopeTransactionExtension ?
# -> http://stackoverflow.com/a/6044925
# -> pyramid_tm (transaction manager) is configured
DBSession = scoped_session(sessionmaker(extension=ZopeTransactionExtension()))

# an unmanaged session
# used by the post() function,
# that needs to commit the session eraly on
DBSession.unmanaged = scoped_session(sessionmaker())


def main(global_config, **settings):
    """
    Return a Pyramid(!) WSGI application.
    """
    engine = create_engine(rnaseqlyze.db_url)
    DBSession.configure(bind=engine)
    DBSession.unmanaged.configure(bind=engine)

    config = Configurator(settings=settings)

    config.scan()

    config.add_route('home', '/')
    config.add_route('analyses', '/analyses')
    config.add_route('analysis', '/analyses/{id}')

    config.add_view(route_name='home', renderer='templates/home.pt')
    config.add_view(route_name='analyses', renderer='templates/create.pt')

    for path in 'less', 'css', 'img', 'js':
        config.add_static_view(path, path)

    return config.make_wsgi_app()


from pyramid.events import subscriber
from pyramid.events import BeforeRender
from pyramid.renderers import get_renderer

@subscriber(BeforeRender)
def add_base_template(event):

    base = get_renderer('templates/base.pt').implementation()

    def path(relative):
        return event['request'].route_path('home') + relative

    event.update({
        'base': base,
        'path': path,
        'version': rnaseqlyze.__version__,
    })
