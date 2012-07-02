#!/usr/bin/python

prefix = "/home/pfeifer/.local"
workdir = "/home/pfeifer/data/rna-seqlyze-workdir"

import site, glob
site.addsitedir(glob.glob(prefix + "/lib/python2*/site-packages")[-1])

from rnaseqlyze.web import wsgi
application = wsgi.get_app(workdir)
