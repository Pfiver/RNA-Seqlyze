"""
A collection of security related functions.

In case a check fails, an exception is raised, otherwise None is returned.
"""

def check_valid_filename(name):
    """
    Assert that the passed name doesn't contain
    any "funny" characters (e.g. ../../../../sensitive.txt)
    """

    max_len = 128
    assert len(name) < max_len, "Filename too long"

    import string
    assert set(srr) < set(string.digits + string.letters + '_'), \
            "Only digits, letters and underscore are allowed in filenames"
