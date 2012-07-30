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

    # make sure to be able to delete files created by webapp
    # as user/group www-data/www-data from the command line
    # (as user/group johndoe/www-data)
    import os
    os.umask(0002)

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
        if False: # production
            e, t = arg, type(arg)
            arg = "%s %s" % (t.__name__, e.args)
        else: # debug
            import traceback
            arg = traceback.format_exc(999)
    Response.__init__(self,
        '%s %s\n\n%s' % (self.code, self.title, arg),
        content_type='text/plain', status='%s %s' % (self.code, self.title))

from pyramid.httpexceptions import WSGIHTTPException
WSGIHTTPException.__init__ = _WHE_init
del WSGIHTTPException.__call__
del WSGIHTTPException.prepare
