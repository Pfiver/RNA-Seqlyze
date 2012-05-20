#!/usr/bin/python

import site
site.addsitedir("/home/pfeifer/.local/lib/python2.6/site-packages")
from trac.web.main import dispatch_request

def application(environ, start_request):
	environ['trac.env_path'] = '/home/pfeifer/data/bt/src/trac_env'
	environ['SCRIPT_NAME'] = environ['SCRIPT_NAME'][:-5]
	return dispatch_request(environ, start_request)
