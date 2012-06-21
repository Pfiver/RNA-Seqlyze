"""
This module implements all kinds of security related checks.

In case a checks fails, an exception is raised, otherwise None is returned.
"""

def check_valid_filename(name):
    """
    Make sure the name doesn't contain
    any "funny" characters - ie. ../../../../sensitive.txt
    """

    max_len = 128
    assert len(name) < max_len, "Filename too long"

    import string
    assert set(srr) < set(string.digits + string.letters + '_'), \
            "Only digits, letters and underscore are allowed in filenames"
