#!/usr/bin/python

webapp_path = '/home/pfeifer/data/bt/src/rna-seqlyze-web/production.ini'

import site
site.addsitedir("/home/pfeifer/.local/lib/python2.6/site-packages")

import pyramid.paster
def application(environ, start_request):
    # strip the ".wsgi" ending
    environ['SCRIPT_NAME'] = environ['SCRIPT_NAME'][:-5]
    return pyramid.paster.get_app(webapp_path, 'main')(environ, start_request)
