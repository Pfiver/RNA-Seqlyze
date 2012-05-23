#!/usr/bin/python

import site
site.addsitedir("/home/pfeifer/.local/lib/python2.6/site-packages")
from pyramid.paster import get_app

def application(environ, start_request):
    environ['SCRIPT_NAME'] = environ['SCRIPT_NAME'][:-5]
    return get_app('/home/pfeifer/data/bt/src/rna-seqlyze-web/development.ini', 'main')(environ, start_request)
