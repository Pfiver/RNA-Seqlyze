import json, urllib2, time
from trac.web.href import Href
from trac.wiki.macros import WikiMacroBase

def get_time(sec):
    localtime = time.localtime(int(sec))
    tz_mins = time.timezone / 60
    tz_hrs = tz_mins / 60 - localtime.tm_isdst
    tz_mins = tz_mins % 60
    tz = "%+.02d:%02d" % (tz_hrs, abs(tz_mins))
    return time.strftime("%Y-%m-%d %H:%M " + tz, localtime)

class BuildbotStatusMacro(WikiMacroBase):
    """BuildbotStatus macro

    depending on the status of the buildbot,
    display a red or a green div with current revision hash
    """
    url = "$URL$"
    revision = "$Rev$"

    def expand_macro(self, formatter, name, text, args):
        get_link = lambda dir_: formatter.req.href("/browser/src/%s" % dir_) + "?rev=master"
        builds = [ json.load(urllib2.urlopen("http://localhost:8010/json/builders/%s/builds/-1" % slave))
                    for slave in "rna-seqlyze", "rna-seqlyze-web" ]
        out = []
        out.append("<ul>")
        for build in builds:
            status = " ".join(build["text"])
            color = ["#F00", "#0F0"]["success" in status]
            out.append(
"""\
<li><a href="%s">%s</a>: <span style="display: inline-block; padding: 4px; border-radius: 4px; background-color: %s">%s</span><br>
 %s<br>
 %s<br>
</li>
""" % (get_link(build["builderName"]), build["builderName"], color, status, build["blame"][0], get_time(build["times"][1])))

        out.append("</ul>")
        return "\n".join(out)
