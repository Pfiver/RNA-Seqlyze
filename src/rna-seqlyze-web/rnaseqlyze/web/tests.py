import unittest
import transaction

from pyramid import testing

from rnaseqlyze.core.orm import DBSession, Analysis

def setup():
    testing.setUp()
    from sqlalchemy import create_engine
    engine = create_engine('sqlite://')
    DBSession.configure(bind=engine)
    from rnaseqlyze.core.orm import Entity, Analysis
    Entity.metadata.create_all(engine)
    with transaction.manager:
        analysis = MyModel(refseq_ns='NC_1234')
        DBSession.add(analysis)

def teardown():
    DBSession.remove()
    testing.tearDown()

def test_it(self):
    from .views import main
    request = testing.DummyRequest()
    info = main(request)
    self.assertEqual(info['analysis'].refseq_ns, 'NC_1234')
