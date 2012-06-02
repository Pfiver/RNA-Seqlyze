#!/usr/bin/python
# encoding: utf-8

import subprocess

from pyramid.view import view_config
from pyramid.config import Configurator
from pyramid.response import Response

child = None

@view_config('start', request_method='GET')
def start(request):
    global child
    child = subprocess.Popen("sleep 10; ls -la", stdout=subprocess.PIPE, shell=True)
    return Response("OK")

@view_config('results', request_method='GET', renderer='json')
def results(request):
    if child:
        res = child.stdout.read()
        child.wait()
        return Response({ "res" : res }, content_type="text/plain")
        return { "res" : res }
    return Response({ "ans" : "start first" }, content_type="text/plain")
    return { "ans" : "start first" }

# scan for @view_config and @subscriber decorators
config = Configurator()
config.scan()

app = config.make_wsgi_app()

#from wsgiref.simple_server import make_server
#server = make_server('0.0.0.0', 8080, app)
#server.serve_forever()

from waitress import serve
serve(app, host='0.0.0.0', port=8080)
