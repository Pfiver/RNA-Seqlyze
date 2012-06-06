# based on http://code.activestate.com/recipes/146306/

def urlopen_multipart(url, data):
	import urllib2
	req = urllib2.Request(url, data)
	try:
		boundary = data[2:data.index("\r")]
	except ValueError, e:
		raise Exception("couldn't find boundary string in data", e)
	req.add_header('Content-Type', 'multipart/form-data; boundary=%s' % boundary)
	return urllib2.urlopen(req)

def urlencode_multipart(fields, files=None):
	"""
	fields is a sequence of (name, value) elements for regular form fields.
	files is a sequence of (name, filename, value) elements for data to be uploaded as files
	Return "multipart/form-data" encoded fields + files
	"""
	import uuid
	boundary = str(uuid.uuid4())
	data = []
	for (key, value) in fields:
		data.append('--' + boundary)
		data.append('Content-Disposition: form-data; name="%s"' % key)
		data.append('')
		data.append(value)
	if files:
		for (key, filename, value) in files:
			data.append('--' + boundary)
			data.append('Content-Disposition: form-data; name="%s"; filename="%s"' % (key, filename))
			data.append('Content-Type: %s' % get_content_type(filename))
			data.append('')
			data.append(value)
	data.append('--' + boundary + '--')
	data.append('')

	return '\r\n'.join(data)

def get_content_type(filename):
	import mimetypes
	return mimetypes.guess_type(filename)[0] or 'application/octet-stream'

if __name__ == "__main__":
	import socket
	import threading
	sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
	sock.bind(('localhost', 9123))
	sock.listen(1)
	def run():
		con = sock.accept()[0]
		print con.recv(4096)
		con.send("200 OK\r\n\r\n")
		con.close()
	threading.Thread(target=run).start()

	data = urlencode_multipart((("123", "456"),))
	urlopen_multipart("http://localhost:9123/asd", data)
