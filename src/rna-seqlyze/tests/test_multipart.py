from nose.tools import *

def test_multipart():
    import io, socket, threading
    def run():
        sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        sock.bind(('localhost', 9123))
        sock.listen(1)
        con = sock.accept()[0]
        con.settimeout(0.1)
        while True:
            try:
                received.write(con.recv(4096).decode("utf-8"))
            except socket.timeout:
                break
        con.send("200 OK\r\n\r\n")
        con.close()
        sock.close()
    received = io.StringIO()
    threading.Thread(target=run).start()

    from rnaseqlyze.multipart \
        import urlencode, urlopen
    data = urlencode((("foo", "bar"),))
    rsp = urlopen("http://localhost:9123/foobar", data)

    received = received.getvalue()
    assert_equal(rsp.getcode(), 200)
    assert_in('POST /foobar HTTP/1.1', received)
    assert_in('Content-Type: multipart/form-data', received)
    assert_in('Content-Disposition: form-data; name="foo"', received)

expected = """\
POST /foobar HTTP/1.1
Accept-Encoding: identity
Content-Length: 133
Host: localhost:9123
Content-Type: multipart/form-data; boundary=d85867df-c194-4e6a-acb4-35bd19ce263f
Connection: close
User-Agent: Python-urllib/2.6

--d85867df-c194-4e6a-acb4-35bd19ce263f
Content-Disposition: form-data; name="foo"

bar
--d85867df-c194-4e6a-acb4-35bd19ce263f--
"""
