# api-doc:
# http://docs.pylonsproject.org/projects/pyramid/en/1.3-branch/narr/hooks.html#registering-configuration-decorators
#
# explanation why it doesn't exist already:
# http://docs.pylonsproject.org/projects/pyramid/en/1.3-branch/designdefense.html#routes-need-relative-ordering
#
# adapted from class view_config(object):
#  https://github.com/Pylons/pyramid/blob/15e3b19/pyramid/view.py#L149
# and groundhog:
#  https://github.com/Pylons/groundhog/blob/master/groundhog.py#L80

import venusian
class route(object):
    def __init__(self, pattern, **view_args):
        self.pattern = pattern
        self.view_args = view_args
    def __call__(self, wrapped):
        venusian.attach(wrapped, self.callback)
        return wrapped
    def callback(self, scanner, name, wrapped):
        scanner.config.add_route(name=wrapped.__name__, pattern=self.pattern)
        scanner.config.add_view(view=wrapped, route_name=wrapped.__name__, **self.view_args)

from wsgiref.simple_server import make_server

from pyramid.view import view_config
from pyramid.config import Configurator
from pyramid.response import Response

@route('/start/{number}', request_method='GET')
def start(request):
    return Response("Number " + request.matchdict["number"])

def main():
    config = Configurator()
    config.scan()
    app = config.make_wsgi_app()
    server = make_server('0.0.0.0', 8080, app)
    server.serve_forever()

main()
