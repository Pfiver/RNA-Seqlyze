def application(environ, start_response):
    status = '200 OK'
    output = 'Hello World!'

    response_headers = [('Content-type', 'text/plain')]
#    response_headers = [('Content-type', 'text/plain'),
#                        ('Content-Length', str(len(output)))]

    start_response(status, response_headers)

#    import subprocess
#    p = subprocess.Popen("ps axfu", stdout=subprocess.PIPE, shell=True)

#    return p.stdout

#    return ("%s: %s\n" % (n,v) for n,v in environ.iteritems())

    return [output]
