#!/usr/bin/python

import site
site.addsitedir("/home/pfeifer/.local/lib/python2.6/site-packages")
from pyramid.paster import get_app

def application(env, start_request):
    pi = len(env['REQUEST_URI']) - len(env['PATH_INFO']) 
    ni = env['REQUEST_URI'].rindex("/",0,pi) + 1
    env['SCRIPT_NAME'] = env['REQUEST_URI'][:pi]
    appname = env['REQUEST_URI'][ni:pi]
    dev_ini = '/home/pfeifer/webapps/'+appname+'/development.ini'
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
    return get_app('/home/pfeifer/webapps/'+appname+'/development.ini', 'main')(env, start_request)
