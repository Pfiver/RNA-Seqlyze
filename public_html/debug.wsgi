#!/usr/bin/python

import os, sys, site, signal, subprocess, trac
from pprint import pformat

site.addsitedir("/home/pfeifer/.local/lib/python2.6/site-packages")
os.environ["PYTHON_EGG_CACHE"] = "/home/pfeifer/.wsgi-egg-cache"

def application(env, start_response):
    def flt(n):
         return not n.startswith("__") and n != "application"
    funcs = filter(flt, globals())
    global environ
    environ = env
    start_response('200 OK', [('Content-type', 'text/plain')])
    if env["QUERY_STRING"] not in funcs:
        nam = os.path.basename(env["SCRIPT_NAME"])
        return [ "usage: %s [%s]" % (nam, "|".join(funcs)) ]
    return globals()[env["QUERY_STRING"]]()

def ps_afux():
    p = subprocess.Popen(["ps", "axfu"], stdout=subprocess.PIPE)
    return p.stdout.readlines()

def sig_int():
    os.kill(os.getpid(), signal.SIGINT)
    return [ "pid: %d\n" % os.getpid() ]

def ps_afux():
    p = subprocess.Popen("ps axfu", stdout=subprocess.PIPE, shell=True)
    return p.stdout

def sys_path():
    return [ "sys.path:\n", pformat(sys.path), "\n"]

def egg_cache():
    return [ os.environ["PYTHON_EGG_CACHE"], "\n" ]

def trac_version():
    return [ trac.__version__, "\n" ]

def printenv():
    return [
        "pid: %d\n" % os.getpid(),
        "environ:\n", pformat(environ), "\n"
    ]
