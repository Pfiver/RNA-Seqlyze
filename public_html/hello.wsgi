def application(environ, start_response):
    status = '200 OK'
    output = 'Hello World!'

    response_headers = [('Content-type', 'text/plain')]
#    response_headers = [('Content-type', 'text/plain'),
#                        ('Content-Length', str(len(output)))]

    start_response(status, response_headers)

#    import subprocess
#    p = subprocess.Popen("ps axfu", stdout=subprocess.PIPE, shell=True)
#
#    return p.stdout

    out = []
    import sys
    from pprint import pformat
    out.append("old path:")
    out.append(pformat(sys.path))

    import site
    site.addsitedir("/home/pfeifer/.local/lib/python2.6/site-packages")

    out.append("new path:")
    out.append(pformat(sys.path))

    import trac
    out.append(str(trac.__version__))

    return [ "\n".join(out) ]

    return ("%s: %s\n" % (n,v) for n,v in environ.iteritems())

    return [output]
