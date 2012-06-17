import logging
log = logging.getLogger(__name__)

from pyramid.config import Configurator
from pyramid.response import Response
from pyramid.httpexceptions import HTTPBadRequest

import rnaseqlyze
from .core import Manager
from .core import AnalysisAlreadyStartedException


def wsgi_app(global_config, **settings):
    """
    Return a Pyramid(!) WSGI application.
    """

    Waitress.manager = Manager()

    config = Configurator(settings=settings, root_factory=Waitress)
    config.add_view(Waitress.text_view, renderer='string')
    return config.make_wsgi_app()


class Waitress(object):

    def __init__(self, request):
        self.analyses = None
        self.analysis = None
        self.status = None
        self.paths = self.pathiterator()

    def __getitem__(self, key):
        self.key = key
        return self.paths.next()

    def pathiterator(self):
        # /analyses
        if self.key != 'analyses':
            raise KeyError()
        self.analyses = True
        yield self

        # /analyses/32
        self.analysis = self.manager.get_analysis(int(self.key))
        if not self.analysis:
            raise KeyError()
        yield self

        # /analyses/32/status
        if self.key != 'status':
            raise KeyError()
        self.status = True
        yield self

        raise KeyError()

    @classmethod
    def text_view(cls, context, request):

        if request.method == 'START':

            if not context.analysis:
                raise HTTPBadRequest("START must be called on /analyses/#")

            try:
                cls.manager.start_analysis(context.analysis)
            except AnalysisAlreadyStartedException, e:
                raise HTTPBadRequest(e)

            return "started analysis #%d" % context.analysis.id

        import pprint
        return pprint.pformat({
            'context': context, # Waitress
            'manager': cls.manager, # Manager
            'analysis': context.analysis, # Analysis
        })


# monkey-patch some HTTPException classes to get simpler error messages

from pyramid.response import Response
def _WHE_init(self, arg=None):
    Exception.__init__(self, arg)
    if isinstance(arg, Exception):
        e, t = arg, type(arg)
        arg = "%s.%s(%s)" % (t.__module__, t.__name__, e.args[0])
    Response.__init__(self, arg + '\n',
        content_type='text/plain', status='%s %s' % (self.code, self.title))

from pyramid.httpexceptions import WSGIHTTPException
WSGIHTTPException.__init__ = _WHE_init
del WSGIHTTPException.__call__
del WSGIHTTPException.prepare
