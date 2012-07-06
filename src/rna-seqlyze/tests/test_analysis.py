import transaction
from sqlalchemy import create_engine
from sqlalchemy.orm import scoped_session, sessionmaker

import rnaseqlyze
from rnaseqlyze.core.entities import Entity, Analysis

Session = scoped_session(sessionmaker())

def setup():
    engine = create_engine('sqlite://')
    Entity.metadata.create_all(engine)
    Session.configure(bind=engine)
    analysis = Analysis(org_accession="NC_1234")
    with transaction.manager:
        Session.add(analysis)

def teardown():
    Session.remove()

from nose.tools import *
def basic_tests():

    a = Session.query(Analysis).get(1)
    assert_false(a.genbankfile_uploaded)
    assert_equal(a.genbankfile_name, "NC_1234.gb")
    assert_equal(a.genbankfile_fa_name, "NC_1234.fa")

    a = Analysis()
    assert_false(a.inputfile_uploaded)

    a = Analysis(org_accession="NC_1234")
    assert_equal(a.org_accession, "NC_1234")
    assert_false(a.genbankfile_uploaded)

    a = Analysis(inputfile_name="SRR123.sra")
    assert_equal(a.inputfile_name, "SRR123.sra")
    assert_true(a.inputfile_uploaded)

    a = Analysis(genbankfile_name="NC_1234.gb")
    assert_true(a.genbankfile_uploaded)
    assert_equals(a.genbankfile_name, "NC_1234.gb")
    assert_equals(a.genbankfile_base_name, "NC_1234")
    assert_equals(a.genbankfile_fa_name, "NC_1234.fa")
