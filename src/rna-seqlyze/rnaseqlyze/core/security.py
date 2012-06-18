def check_accession(acc):
    """
    Make sure the accession doesn't contain
    any "funny" characters - ie. ../../../../sensitive.txt
    """

    assert len(acc) < 50

    valid_chars = [chr(c) for c in
            range(ord("A"), ord("Z") + 1)
            + range(ord("0"), ord("9") + 1)] + ['_', '.']

    for c in acc:
        assert c in valid_chars, "Bad character in accession number: '%c'" % c
