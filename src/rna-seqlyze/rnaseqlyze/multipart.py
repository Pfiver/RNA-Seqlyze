"""
Multipart form-data handling
based on http://code.activestate.com/recipes/146306/
"""
import uuid, urllib2, mimetypes

def urlopen(url, data=None):
    if isinstance(url, basestring):
        rq = urllib2.Request(url, data)
    elif isinstance(url, urllib2.Request):
        rq = url
        data = rq.data
    else:
        raise Exception("'url' parameter must be a string or urllib2.Request")

    try:
        boundary = data[2:data.index("\r")]
    except ValueError, e:
        raise Exception("couldn't find boundary string in data", e)
    rq.add_header('Content-Type', 'multipart/form-data; boundary=%s' % boundary)
    return urllib2.urlopen(rq)

def urlencode(fields, files=None):
    """
    :param asd:
      is a sequence of ``(name, value)`` elements for regular form fields.

    :param files:
      is a sequence of ``(name, filename, value)``
      elements for data to be uploaded as files

    :returns:
      ``str`` of **multipart/form-data** encoded fields + files
    """
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
            data.append('Content-Disposition: form-data' \
                            '; name="%s"; filename="%s"' % (key, filename))
            data.append('Content-Type: %s' % get_content_type(filename))
            data.append('')
            data.append(value)
    data.append('--' + boundary + '--')
    data.append('')

    return '\r\n'.join(data)

def get_content_type(filename):
    return mimetypes.guess_type(filename)[0] or 'application/octet-stream'
