#!/usr/bin/python

import site
site.addsitedir("/home/pfeifer/.local/lib/python2.6/site-packages")

def application(env, start_response):
    def flt(n):
         return not n.startswith("__") and n != "application"
    funcs = filter(flt, globals())
    global environ
    environ = env
    start_response('200 OK', [('Content-type', 'text/plain')])
    if env["QUERY_STRING"] not in funcs:
        import os
        nam = os.path.basename(env["SCRIPT_NAME"])
        return [ "usage: %s [%s]" % (nam, "|".join(funcs)) ]
    return globals()[env["QUERY_STRING"]]()

def sig_int():
    import os, signal
    pid = os.getpid()
    os.kill(pid, signal.SIGINT)
    return [ "signaled pid %d\n" % pid ]

def ps_afux():
    import subprocess
    p = subprocess.Popen(["ps", "axfu"], stdout=subprocess.PIPE)
    return p.stdout

def sys_path():
    import pprint, sys
    return [ "sys.path:\n", pprint.pformat(sys.path), "\n"]

def trac_version():
    import trac
    return [ trac.__version__, "\n" ]

def printenv():
    import pprint, os
    return [
        "pid: %d\n" % os.getpid(),
        "environ:\n", pprint.pformat(environ), "\n",
        "os.environ:\n", pprint.pformat(dict(os.environ)), "\n",
    ]
