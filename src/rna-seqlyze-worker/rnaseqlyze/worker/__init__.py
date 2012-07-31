"""
**pyramid.worker** is a Pyramid Web Framework Application.

The framework is used here to keep things simple. Even thought not many of the
frameworks features are used, building this "-worker" part of the project as a
Pyramid Web Framework Application, just like the "-web" part, hopefullly makes
it easy to understand for anybody already understanding the "-web" part.

The key features from then Pyramid Web Framework used here, are

 1) The "pserve" command, which makes running the application as a unix daemon
    process very simple. It's direct use has actually been depreciated during
    the development and a custom command, "rnas-worker", has been created,
    which uses the same python functions and modules like "pserve".

 2) The pyramid.config.Configurator class, that is used to define the
    applications "routes" and "view callables". These "views" provide the
    applications interface.  They are served on a tcp port bound to localhost
    (127.0.0.1) and are therefore only available to processes running on the
    same host.
    
The "-worker" applications interface has "HTTP-like" semantics.

The following commands are accepted:

 - ``GET /analyses/{id}``:      Show the current status.
 - ``START /analyses/{id}``:    Start processing an analysis.

Only for development purposes, one additional command exists:

 - ``RESTART /analyses/{id}``:  Restart an analysis
                                that has already been started.

The commands are executed by the "-web" part of the application by subclassing
HTTPRequest and coverriding the get_method() function. They can also be executed
from the command line however, using the popular "curl" binary with the "-X"
option, e.g.

 - ``curl -X GET localhost:/analyses/3``
 - ``curl -X START localhost:/analyses/3``
 - ``curl -X RESTART localhost:/analyses/3``

"""

import logging
log = logging.getLogger(__name__)

from pyramid.view import view_config
from pyramid.view import view_defaults
from pyramid.config import Configurator
from pyramid.response import Response
from pyramid.httpexceptions import (
        HTTPError,
        HTTPBadRequest,
        HTTPInternalServerError,
)

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker, scoped_session
from zope.sqlalchemy import ZopeTransactionExtension

import rnaseqlyze
from rnaseqlyze.core.entities import Analysis
from rnaseqlyze.worker.core import (
        Manager,
        ManagerBusyException,
        AnalysisAlreadyStartedException,
)

DBSession = scoped_session(sessionmaker(extension=ZopeTransactionExtension()))

def main(global_config, **settings):
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
    def status(self):
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
def error_view(request):
    errdict = {
        AnalysisAlreadyStartedException:    HTTPRNASBadRequest,
        ManagerBusyException:               HTTPRNASWorkerError,
    }
    type_ = type(request.exc_info[1])
    return errdict.get(type_, HTTPRNASWorkerError)(request.exc_info)

class HTTRNASError(HTTPError):
    def __call__(self, environ, start_response):
        # steer clear from WSGIHTTPException.__call__,
        # which calls WSGIHTTPException.prepare,
        # which sets content_type = 'text/html'
        return Response.__call__(self, environ, start_response)
    def __init__(self, error):
        err = exc_info[1]
        log.error(repr(err))
        cls = e.__class__.__name__
        if not e.args:
            explanation = "%s" % cls
        else:
            explanation = "%s: %s" % (cls, e.args[0])
        import traceback
        detail = (str(err) +
                  '\n\nStack trace:\n' +
                  ''.join(traceback.format_tb(exc_info[2])))
        log.debug(detail)
        HTTPError.__init__(self, detail, content_type='text/plain')

class HTTRNASBadRequest(HTTRNASError):
    code = 400
    title = "Bad Request"

class HTTRNASWorkerError(HTTRNASError):
    code = 500
    title = "RNA-Seqlyze Worker Error"
