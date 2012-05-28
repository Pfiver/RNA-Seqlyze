#!/usr/bin/python

trac_env = '/home/pfeifer/data/bt/src/trac_env'

import site
site.addsitedir("/home/pfeifer/.local/lib/python2.6/site-packages")

import trac.web.main
def application(environ, start_request):
    # strip the ".wsgi" ending
    environ['SCRIPT_NAME'] = environ['SCRIPT_NAME'][:-5]
    environ['trac.env_path'] = trac_env
    return trac.web.main.dispatch_request(environ, start_request)
