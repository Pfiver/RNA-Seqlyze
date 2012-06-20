# based on http://code.activestate.com/recipes/146306/

def urlopen(url, data):
	import urllib2
	req = urllib2.Request(url, data)
	try:
		boundary = data[2:data.index("\r")]
	except ValueError, e:
		raise Exception("couldn't find boundary string in data", e)
	req.add_header('Content-Type', 'multipart/form-data; boundary=%s' % boundary)
	return urllib2.urlopen(req)

def urlencode(fields, files=None):
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
