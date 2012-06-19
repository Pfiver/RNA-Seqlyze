import logging
log = logging.getLogger(__name__)

from pyramid.view import view_config
from pyramid.view import view_defaults
from pyramid.config import Configurator
from pyramid.response import Response
from pyramid.httpexceptions import (
        HTTPError, HTTPBadRequest, HTTPInternalServerError
)
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker, scoped_session
from zope.sqlalchemy import ZopeTransactionExtension

import rnaseqlyze
from rnaseqlyze.core.orm import Analysis
from .core import Manager
from .core import AnalysisAlreadyStartedException, ManagerBusyException

DBSession = scoped_session(sessionmaker(extension=ZopeTransactionExtension()))

def wsgi_app(global_config, **settings):
    """
    Return a Pyramid(!) WSGI application.
    """

    engine = create_engine(rnaseqlyze.db_url)
    DBSession.configure(bind=engine)

    Waitress.manager = Manager()

    config = Configurator(settings=settings)
    config.add_route('analyses', '/analyses/{id}')
    config.scan()
    return config.make_wsgi_app()


@view_defaults(route_name='analyses', renderer='string')
class Waitress(object):

    def __init__(self, request):
        id = int(request.matchdict['id'])
        self.analysis = DBSession.query(Analysis).get(id)

    @view_config(request_method='GET')
    def get(self):
        import pprint
        return pprint.pformat({
            'context': self, # Waitress
            'manager': self.manager, # Manager
            'analysis': self.analysis, # Analysis
        })

    @view_config(request_method='START')
    def start(self):
        self.manager.analysis_requested(self.analysis)
        return "started analysis #%d" % self.analysis.id

    @view_config(request_method='RESTART')
    def restart(self):
        self.manager.analysis_requested(self.analysis, True)
        return "restarted analysis #%d" % self.analysis.id


@view_config(context=Exception)
def error_view(error, request):
    errdict = {
        AnalysisAlreadyStartedException:    HTTPBadRequest,
        ManagerBusyException:               HTTPInternalServerError,
    }
    if isinstance(error, HTTPError):
        return error
    elif type(error) in errdict:
        return errdict[type(error)](error)
    else:
        return HTTPInternalServerError(error)


# monkey-patch some HTTPException classes to get simpler error messages

from pyramid.response import Response
def _WHE_init(self, arg=None):
    Exception.__init__(self, arg)
    if isinstance(arg, Exception):
        e, t = arg, type(arg)
        arg = "%s.%s %s" % (t.__module__, t.__name__, e.args)
        import traceback
        arg += '\n' + traceback.format_exc(999)
    Response.__init__(self,
        'HTTP %s %s: %s\n' % (self.code, self.title, arg),
        content_type='text/plain', status='%s %s' % (self.code, self.title))

from pyramid.httpexceptions import WSGIHTTPException
WSGIHTTPException.__init__ = _WHE_init
del WSGIHTTPException.__call__
del WSGIHTTPException.prepare
