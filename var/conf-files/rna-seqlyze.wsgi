#!/usr/bin/python

workdir = "@@WORKDIR@@"

from rnaseqlyze.web import wsgi
application = wsgi.get_app(workdir)
