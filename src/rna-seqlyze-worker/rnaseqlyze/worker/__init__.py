import logging
log = logging.getLogger(__name__)

from pyramid.config import Configurator
from pyramid.response import Response

import rnaseqlyze
from .core import Manager


def wsgi_app(global_config, **settings):
    """
    Return a Pyramid(!) WSGI application.
    """

    Waitress.manager = Manager()

    config = Configurator(settings=settings, root_factory=Waitress)
    config.add_view(Waitress.json_view, renderer='json')
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
    def json_view(cls, context, request):
        if request.method == 'START':
            if context.analysis:
                n = cls.manager.start_analysis(context.analysis)
            elif context.analyses:
                n = cls.manager.start_all_analyses()
            else:
                raise Exception("start what?")
            return Response("started %d analyses" % n)

        return {
            'context': str(context), # Waitress
            'manager': str(cls.manager), # Manager
            'analysis': str(context.analysis), # Analysis
#            'started': context.analysis.started,
        }
