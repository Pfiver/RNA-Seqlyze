#!/usr/bin/python

trac_env = '@@TOPDIR@@/var/trac-env'

import site
site.addsitedir("@@PREFIX@@/lib/python2.6/site-packages")

import os
os.umask(0002)

import trac.web.main
def application(environ, start_request):
    # strip the ".wsgi" ending
    environ['SCRIPT_NAME'] = environ['SCRIPT_NAME'][:-5]
    environ['trac.env_path'] = trac_env
    return trac.web.main.dispatch_request(environ, start_request)
