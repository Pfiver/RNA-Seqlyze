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
    assert len(name) < max_len, \
        "Filename lengh (%d) exceeds limit of %d" % (len(name), max_len)

    valid_chars = [chr(c) for c in
            range(ord("A"), ord("Z") + 1)
            + range(ord("a"), ord("z") + 1)
            + range(ord("0"), ord("9") + 1)] + ['_', '.']

    for c in name:
        assert c in valid_chars, "Bad character in filename: '%c'" % c
