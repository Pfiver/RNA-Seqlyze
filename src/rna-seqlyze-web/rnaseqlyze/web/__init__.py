"""
**pyramid.web** is a Pyramid Web Framework Application.

To learn more about the applications architecture, head over to the wonderful
world of the pyramid web framework at http://www.pylonsproject.org/.

This application has been created using the ``pcreate`` command with the ``-s
alchemy`` option to create and sqlalchemy scaffold. There is plenty of `very
good documentation <http://docs.pylonsproject.org/projects/pyramid/en/latest/\
narr/project.html#scaffolds-included-with-pyramid>`_ available on how to do it.

In case you have trouble with anything pyramid-related, use the `source code on
github <https://github.com/Pylons/pyramid>`_ or ask 'mcdonc' on freenode irc
`#pyramid <http://webchat.freenode.net/?channels=#pyramid>`_.
"""

import logging
log = logging.getLogger(__name__)

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker, scoped_session
from pyramid.config import Configurator
#: the zope transaction extension
from zope.sqlalchemy import ZopeTransactionExtension

import rnaseqlyze
from rnaseqlyze.web.jsonx import jsonx

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
    Create and return a Pyramid WSGI application.
    """
    log.debug("rnaseqlyze.web version %s : main()" % rnaseqlyze.__version__)

    # make sure to be able to delete files created by webapp
    # as user/group www-data/www-data from the command line
    # (as user/group johndoe/www-data)
    import os
    os.umask(0002)

    engine = create_engine(rnaseqlyze.db_url)
    DBSession.configure(bind=engine)
    DBSession_unmanaged.configure(bind=engine)

    config = Configurator(settings=settings)

    config.add_renderer('jsonx', jsonx)

    config.scan()

    config.add_route('home', '/')
    config.add_route('upload', '/upload')
    config.add_route('analyses', '/analyses')
    config.add_route('analysis', '/analyses/{id}')
    config.add_route('analysis_files', '/analyses/{id}/files*subpath')

    config.add_route('analysis_rest', '/rest/analyses/{id}')
    config.add_route('analysis_logs_rest', '/rest/analyses/{id}/logs')
    config.add_route('analysis_files_rest', '/rest/analyses/{id}/files')

    config.add_route('organisms_rest', '/rest/organisms')

    for path in 'less', 'css', 'img', 'js':
        config.add_static_view(path, path)

    return config.make_wsgi_app()

from pyramid.events import subscriber
from pyramid.events import BeforeRender
from pyramid.renderers import get_renderer

@subscriber(BeforeRender)
def before_render(event):
    """
    This function is called by Pyramid after the view callable has returned
    and before the renderer (json / chameleon) is called. We inject some
    convienience functions and objects that are used in the `zope page
    templates <http://pagetemplates.org/docs/latest/reference.html>`_
    (.pt files) which `the chameleon template engine
    <http://pagetemplates.org/>`_ then renders.
    """

    base = get_renderer('templates/base.pt').implementation()

    rq = event['request']
    path = lambda sub: rq.route_path('home') + sub
    relpath = lambda sub: rq.current_route_path() + '/' + sub

    event.update({
        'base': base,
        'path': path,
        'relpath': relpath,
        'version': rnaseqlyze.__version__,
        'debug': log.getEffectiveLevel() <= logging.DEBUG,
    })
