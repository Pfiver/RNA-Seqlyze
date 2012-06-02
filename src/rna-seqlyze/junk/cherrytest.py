#!/usr/bin/python
# encoding: utf-8

import cherrypy
from cherrypy import expose

class Converter:
    @expose
    def index(self):
        return "Hello World!"

    @expose
    def fahr_to_celc(self, *vals):
        temp = (float(vals[0]) - 32) * 5 / 9
        return "%.01f" % temp

    @expose
    def celc_to_fahr(self, degrees):
        temp = float(degrees) * 9 / 5 + 32
        return "%.01f" % temp


#from cherrypy.process.plugins import Daemonizer
#d = Daemonizer(cherrypy.engine)
#d.subscribe()

cherrypy.quickstart(Converter())
