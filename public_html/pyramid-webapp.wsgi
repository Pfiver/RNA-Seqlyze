#!/usr/bin/python

webapps_env = '/home/pfeifer/webapps-virtualenv'

import site
site.addsitedir(webapps_env + "/lib/python2.6/site-packages")
site.addsitedir("/home/pfeifer/.local/lib/python2.6/site-packages")

from pyramid.paster import get_app
from paste.script.util import logging_config
def application(env, start_request):

    base_path_len = env['SCRIPT_NAME'].rindex("/") + 1
    try:
        script_name_len = env['REQUEST_URI'].index("/", base_path_len)
    except ValueError:
        script_name_len = len(env['REQUEST_URI'])

    env['SCRIPT_NAME'] = env['REQUEST_URI'][:script_name_len]
    appname = env['SCRIPT_NAME'][base_path_len:]

    basedir = webapps_env + '/apps/' + appname
    dev_ini = basedir + '/development.ini'

    import os
    if not os.path.exists(dev_ini):
        start_request('404 Not Found', [("Content-Type", "text/html")])
        return [
"""
<!DOCTYPE HTML PUBLIC "-//IETF//DTD HTML 2.0//EN">
<html><head>
<title>404 Not Found</title>
</head><body>
<h1>Not Found</h1>
<p>The requested URL %s was not found on this server.</p>
</body></html>
"""
                % env['REQUEST_URI'] ]

    logging_config.fileConfig(dev_ini, {'here': basedir})

    return get_app(dev_ini, 'main')(env, start_request)
