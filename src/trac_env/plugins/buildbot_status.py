import json, urllib2, time
from trac.web.href import Href
from trac.wiki.macros import WikiMacroBase

class BuildbotStatusMacro(WikiMacroBase):
    """BuildbotStatus macro

    depending on the status of the buildbot,
    display a red or a green div with current revision hash
    """
    url = "$URL$"
    revision = "$Rev$"

    def expand_macro(self, formatter, name, text, args):
        build = json.load(urllib2.urlopen("http://localhost:8010/json/builders/runtests/builds/-1"))

        status = build["text"][0]
        color = ["red", "green"][status == "success"]
        return template % (color, status, formatter.req.href("/log"), build["sourceStamp"]["changes"][0]["at"])

template = """\
<ul>
<li>Build Status: <span style="padding: 2px; background-color: %s">%s</span></li>
<li>Latest Changes:<br/><a href="%s">%s</a></li>
</ul>
"""
