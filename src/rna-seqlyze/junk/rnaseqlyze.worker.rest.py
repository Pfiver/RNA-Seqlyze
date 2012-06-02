#!/usr/bin/python

"""
http://lucidcontraptions.com/2010/08/02/rest-service-with-twisted-python/
"""

import sys, json
import subprocess

from twisted.internet import reactor
from twisted.python import log
from twisted.web import server, resource

child = None

class RESTResource(resource.Resource):
    def getChild(self, path, request):
        return NotFound()

class NotFound(RESTResource):
    def render_GET(self, request):
        request.setResponseCode(404)
        return "Not Found\n"

class StartProcessing(RESTResource):
    def render_GET(self, request):
        log.msg("StartProcessing")
        global child
        child = subprocess.Popen("ls -la /", stdout=subprocess.PIPE, shell=True)
        return "OK\n"

class GetResults(RESTResource):
    def render_GET(self, request):
        log.msg("GetResults")
        if child:
            return child.stdout.read()
        return "No results\n"

root = NotFound()
root.putChild("start", StartProcessing())
root.putChild("results", GetResults())

log.startLogging(sys.stdout)
log.msg("Starting server")

service = server.Site(root)
