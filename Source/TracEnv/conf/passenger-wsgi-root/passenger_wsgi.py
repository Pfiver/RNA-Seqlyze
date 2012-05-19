#!/home/aeshoh6c/bin/python2.7
# -*- coding: utf-8 -*-
#
# Copyright (C)2008-2009 Edgewall Software
# Copyright (C) 2008 Noah Kantrowitz <noah@coderanger.net>
# All rights reserved.
#
# This software is licensed as described in the file COPYING, which
# you should have received as part of this distribution. The terms
# are also available at http://trac.edgewall.org/wiki/TracLicense.
#
# This software consists of voluntary contributions made by many
# individuals. For the exact contribution history, see the revision
# history and logs, available at http://trac.edgewall.org/log/.
#
# Author: Noah Kantrowitz <noah@coderanger.net>
import os, urllib

def application(environ, start_request):
    uri = "/"
    path = urllib.unquote(environ["PATH_INFO"][len(environ["REDIRECT_RRDIR"]):])
    if 'HTTPS' in environ:
	uri = path[:path.index("/", 1) + 1]
	path = path[path.index("/", 1):]
#    start_request('200 OK', [('Content-type', 'text/plain')])
#    return [str(
#        uri + "\n" + path + "\n" +
#	"".join("%s -> %s\n" % (var, val) for var, val in environ.iteritems()) +
#	"---------\n" +
#	"".join("%s -> %s\n" % (var, val) for var, val in os.environ.iteritems())
#    )]
    environ["SCRIPT_NAME"] = uri
    environ["PATH_INFO"] = path
    if not 'trac.env_parent_dir' in environ:
        environ.setdefault('trac.env_path', '/home/aeshoh6c/data/pfeifer.ch/RNA-seqlyze')
    if 'PYTHON_EGG_CACHE' in environ:                                           
        os.environ['PYTHON_EGG_CACHE'] = environ['PYTHON_EGG_CACHE']
    elif 'trac.env_path' in environ:
        os.environ['PYTHON_EGG_CACHE'] = \
            os.path.join(environ['trac.env_path'], '.egg-cache')
    elif 'trac.env_parent_dir' in environ:
        os.environ['PYTHON_EGG_CACHE'] = \
            os.path.join(environ['trac.env_parent_dir'], '.egg-cache')
    from trac.web.main import dispatch_request
    return dispatch_request(environ, start_request)
