#!/usr/bin/python

def application(env, start_response):
    def flt(n):
         return not n.startswith("__") and n != "application"
    funcs = filter(flt, globals())
    global environ
    environ = env
    start_response('200 OK', [('Content-type', 'text/plain')])
    if env["QUERY_STRING"] not in funcs:
        import os.path
        nam = os.path.basename(env["SCRIPT_NAME"])
        return [ "usage: %s [%s]" % (nam, "|".join(funcs)) ]
    return globals()[env["QUERY_STRING"]]()

def sig_int():
    import os, signal, subprocess
    p = subprocess.Popen("ps axfu", stdout=subprocess.PIPE, shell=True)
    out = p.stdout.readlines()
    os.kill(os.getpid(), signal.SIGINT)
    p = subprocess.Popen("ps axfu", stdout=subprocess.PIPE, shell=True)
    out += p.stdout.readlines()
    return out

def ps_afux():
    import subprocess
    p = subprocess.Popen("ps axfu", stdout=subprocess.PIPE, shell=True)
    return p.stdout

def sys_path():
    out = []
    import sys, site
    from pprint import pformat
    out.append("old path:")
    out.append(pformat(sys.path))
    site.addsitedir("/home/pfeifer/.local/lib/python2.6/site-packages")
    out.append("new path:")
    out.append(pformat(sys.path))
    return [ "\n".join(out) ]

def trac_version():
    out = []
    import trac
    out.append(str(trac.__version__))
    return [ "\n".join(out) ]

def printenv():
    from pprint import pformat
    return [ pformat(environ) ]
