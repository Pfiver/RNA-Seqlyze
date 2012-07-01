#!/usr/bin/python

# the next two lines can be removed
# if the application is installed under /usr/lib/python2.x
import site
site.addsitedir("/home/pfeifer/.local/lib/python2.6/site-packages")

# when deploying the application, adjust the workdir path here
from rnaseqlyze.web.wsgi import get_app
application = get_app(workdir="/home/pfeifer/data/rna-seqlyze-workdir")
