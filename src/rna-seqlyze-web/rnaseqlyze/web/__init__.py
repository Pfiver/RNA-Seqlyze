"""
**pyramid.web** is a pyramid web framework application.

It has initially been created using the ``pcreate`` command
with the ``-s alchemy`` option to create and sqlalchemy
scaffold. There is plenty of `very good documentation
<http://docs.pylonsproject.org/projects/pyramid/en/latest/"""
"""narr/project.html#scaffolds-included-with-pyramid>`_
available on how to do it.
"""

import logging
log = logging.getLogger(__name__)

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker, scoped_session
from pyramid.config import Configurator
#: the zope transaction extension
from zope.sqlalchemy import ZopeTransactionExtension

import rnaseqlyze

#: a session managed by
#: ZopeTransactionExtension
#:
#: - http://stackoverflow.com/a/6044925
#: - pyramid_tm (transaction manager) is configured
DBSession = scoped_session(sessionmaker(extension=ZopeTransactionExtension()))

#: an unmanaged session
#:
#: used by :meth:`rnaseqlyze.web.views.post`
#: because the session needs to be commited early there
DBSession_unmanaged = scoped_session(sessionmaker())

def main(global_config, **settings):
    """
    Return a Pyramid(!) WSGI application.
    """
    # make sure to be able to delete files created by webapp
    # as user/group www-data/www-data from the command line
    # (as user/group johndoe/www-data)
    import os
    os.umask(0002)

    engine = create_engine(rnaseqlyze.db_url)
    DBSession.configure(bind=engine)
    DBSession_unmanaged.configure(bind=engine)

    config = Configurator(settings=settings)

    config.scan()

    config.add_route('home', '/')
    config.add_route('upload', '/upload')
    config.add_route('analyses', '/analyses')
    config.add_route('analysis', '/analyses/{id}')
    config.add_route('analysis_files', '/analyses/{id}/files*subpath')

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
