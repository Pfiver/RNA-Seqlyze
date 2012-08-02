import cgi, threading, BaseHTTPServer
class Server(threading.Thread):
    class Handler(BaseHTTPServer.BaseHTTPRequestHandler):
        def do_POST(self):
            form = cgi.FieldStorage(
                fp=self.rfile,
                headers=self.headers,
                environ={'REQUEST_METHOD': 'POST'})
            for key in form:
                print >> self.wfile, key, form[key].value,
    def __init__(self):
        threading.Thread.__init__(self)
        self.hts = BaseHTTPServer.HTTPServer(
                    ('localhost', 9123), self.Handler)
    def run(self):
        self.hts.handle_request()

from nose.tools import *
def test_multipart():
    Server().start()
    from rnaseqlyze.multipart \
        import urlencode, urlopen
    data = urlencode((("foo", "bar"),))
    rsp = urlopen("http://localhost:9123/", data)
    assert_equals(rsp.getcode(), 200, "bad response code")
    assert_equals(rsp.read(), "foo bar", "bad response from server")
