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
    if child:
        return Response("In Progress")
    else:
        child = subprocess.Popen("sleep 10; ls -la", stdout=subprocess.PIPE, shell=True)
        return Response("Started")

@view_config('results', request_method='GET', renderer='string')
def results(request):
    global child
    if child:
        stdout = child.stdout
        if child.poll() != None:
            child.wait()
            child = None
        return Response(body_file=stdout, content_type='text/plain')
    return "start first\n"

def main(global_config, **settings):
    config = Configurator(settings=settings)
    config.scan()
    return config.make_wsgi_app()
