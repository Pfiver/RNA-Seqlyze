"""
Pyramid JSON renderer that serializes arbitrary objects

Copies the object's __dict__, looks up all attrs
in all base classes's __dict__'s on the object
and then strips any unknown attribute types.
"""

#import logging
#log = logging.getLogger(__name__)

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
    attrs = dict((attr, getattr(obj, attr))
                    for base in obj.__class__.__bases__
                    for attr in base.__dict__
                    if attr[0] != '_')

    attrs.update(obj.__dict__)
#    log.debug(attrs)
    return dict(filter(filter_attributes, attrs.iteritems()))

none_type = type(None)
def filter_attributes(kv):
    """
    Helper function for render_object
    """
    if kv[0][0] != '_' and \
       type(kv[1]) in (none_type, bool, int, long, float, str, unicode, list):
        return True
    return False
