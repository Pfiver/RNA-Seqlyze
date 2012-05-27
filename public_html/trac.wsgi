#!/usr/bin/python

import os, site
site.addsitedir("/home/pfeifer/.local/lib/python2.6/site-packages")
os.environ["PYTHON_EGG_CACHE"] = "/home/pfeifer/.wsgi-egg-cache"

from trac.web.main import dispatch_request
def application(environ, start_request):
    environ['SCRIPT_NAME'] = environ['SCRIPT_NAME'][:-5]
    environ['trac.env_path'] = '/home/pfeifer/data/bt/src/trac_env'
    return dispatch_request(environ, start_request)
