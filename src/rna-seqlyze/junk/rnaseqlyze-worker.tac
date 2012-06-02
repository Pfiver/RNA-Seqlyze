#!/usr/bin/python
# encoding: utf-8

"""
This is the rna-seqlyze unix worker daemon, which does
all the heavy-lifting on behalf of the web fronteb.

It is built around the "twisted" networking engine.
http://twistedmatrix.com/

To run the daemon on a unix system is has to be
started by SysV-init / upstart / systemd, or whatever
system the undelying distribution is using supports. The
command line to start a daemon that detaches itself from the
console, would be

    $ twistd -y rnaseqlyze-worker.tac

Once started, the daemon listens on port 127.0.0.0:5877
for commands from a local client, i.e. rna-seqlyze-web.

"""

import os

from twisted.application import service, internet
from twisted.web import static, server

from rnaseqlyze.worker import rest

# This object -- 'application' -- is what twistd will look for.
application = service.Application("rna-seqlyze-worker")

# create the sole service this application provides
service = internet.TCPServer(8080, rest.service)

# attach the service to its the application
service.setServiceParent(application)
