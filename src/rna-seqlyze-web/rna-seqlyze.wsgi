#!/usr/bin/python

# Usage:
#       - configure the two directories below
#       - move this file to /var/www/.../rna-seqlyze.wsgi
#       - configure apache to execute the file with mod_wsgi

prefix = "/home/biopython/.local"
workdir = "/home/biopython/data/rna-seqlyze-workdir"

import site, glob
site.addsitedir(glob.glob(prefix + "/lib/python2*/site-packages")[-1])

from rnaseqlyze.web import wsgi
application = wsgi.get_app(workdir)
