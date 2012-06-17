import transaction
from pyramid import testing

from rnaseqlyze.web import DBSession
from rnaseqlyze.core.orm import Analysis

def setup():
    testing.setUp()
    from sqlalchemy import create_engine
    engine = create_engine('sqlite://')
    DBSession.configure(bind=engine)
    from rnaseqlyze.core.orm import Entity, Analysis
    Entity.metadata.create_all(engine)
    with transaction.manager:
        analysis = Analysis()
        DBSession.add(analysis)

def teardown():
    DBSession.remove()
    testing.tearDown()

from nose.tools import *
def most_basic_test():
    from rnaseqlyze.web.views import display
    request = testing.DummyRequest()
    request.matchdict['id'] = "1"
    info = display(request)
    assert_equal(info['analysis'].id, 1)
