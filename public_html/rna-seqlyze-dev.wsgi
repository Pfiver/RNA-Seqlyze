#!/usr/bin/python

webapp_path = '/home/pfeifer/data/bt/src/rna-seqlyze-web/'

import os, site
site.addsitedir("/home/pfeifer/.local/lib/python2.6/site-packages")
os.environ["PYTHON_EGG_CACHE"] = "/home/pfeifer/.wsgi-egg-cache"

from pyramid.paster import get_app
def application(environ, start_request):
    # strip the ".wsgi" ending
    environ['SCRIPT_NAME'] = environ['SCRIPT_NAME'][:-5]
    return get_app(webapp_path + 'development.ini', 'main')(environ, start_request)
