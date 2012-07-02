"""
Pyramid JSON renderer that serializes arbitrary objects
by copying the __dict__ and stripping any unknown attributes
"""

import json

#: a custom json renderer
jsonx = lambda info: render_json

def render_json(value, system):
    """
    custom json renderer implementation

    based on http://git.io/a6BFGQ#L169
    """

    request = system.get('request')
    if request is not None:
        response = request.response
        response.content_type = 'application/json'
    return json.dumps(value, default=render_object, indent=4)

def render_object(obj):
    """
    "default" function for json.dumps()
    """
    return dict(filter(filter_attributes, obj.__dict__.iteritems()))

none_type = type(None)
def filter_attributes(kv):
    """
    Helper function for render_object
    """
    if kv[0][0] != '_' and \
       type(kv[1]) in (none_type, bool, int, long, float, str, unicode, list):
        return True
    return False
