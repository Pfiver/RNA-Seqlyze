#!/usr/bin/python
# encoding: utf-8

from twisted.internet import protocol, reactor, ssl
from twisted.protocols import amp
from twisted.python import log

import sys

cert = ssl.PrivateCertificate.loadPEM(open('ssl.pem').read())
certParams = {'tls_localCertificate': cert}

class SupCommand(amp.Command):
    pass

class ServerProtocol(amp.AMP):
    @SupCommand.responder
    def sup(self):
        log.msg('got sup')
        return {}

    @amp.StartTLS.responder
    def startTLS(self):
        return certParams

class ServerFactory(protocol.ServerFactory):
    protocol = ServerProtocol

class ClientProtocol(amp.AMP):
    def connectionMade(self):
        d = self.callRemote(amp.StartTLS, **certParams)
        d.addCallback(lambda ign: self.callRemote(SupCommand))

class ClientFactory(protocol.ClientFactory):
    protocol = ClientProtocol

if sys.argv[1] == 'client':
    reactor.connectTCP('localhost', 9991, ClientFactory())
elif sys.argv[1] == 'server':
    reactor.listenTCP(9991, ServerFactory())

log.startLogging(sys.stdout)
reactor.run()
