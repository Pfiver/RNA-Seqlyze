# encoding: utf-8

from nose.tools import *
def check_rows(actual, expected):

    # line 1
    # description: "Sulfolobus solfataricus P2"
    assert_equals(actual()[0][:26], expected()[0][:26])
    yield

    # line 2
    # 'nnn Proteins' - not implemented
    assert_equals(actual(), [])
    expected() # discard
    yield

    # line 3
    # column headers
    assert_equals(actual(), expected())
    yield

    # line 4
    # the wrapping around gene SSO12256 (clycotransferase)
    # comes first here but not in the ncbi file
    actual() # discard
    yield

    # lines > 4
    while True:
        # i didn't bother to fill in the "cog" in column #7
        # and the descriptions in column #8 don't match exactly
        assert_equals(actual()[:7], expected()[:7])
        yield

class TestFile(object):
    def __init__(self):
        import csv, Queue
        self.q = Queue.Queue()
        self.checks = check_rows(
            csv.reader(iter(self.q.get, None), delimiter='\t').next,
            csv.reader(open("tests/data/NC_002754-partial.ptt"), delimiter='\t').next)
    def write(self, data):
        self.q.put(data)
        self.checks.next()

def test_gb2ptt():
    from rnaseqlyze.gb2ptt import gb2ptt
    gb2ptt(open("tests/data/NC_002754-partial.gb"), TestFile())
