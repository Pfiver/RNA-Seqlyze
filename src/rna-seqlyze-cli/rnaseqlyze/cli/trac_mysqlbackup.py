def main(argv=sys.argv):

    # http://docs.python.org/library/tempfile.html?highlight=tempfile#tempfile.mkstemp

    from tempfile import mkstemp
    tmp_fd, tmp_path = mkstemp()

    import os
    tmp_file = os.fdopen(tmp_fd, "w")

    # trac-admin . config get trac database |
    #     python -c 'import sys, urlparse; parts=urlparse.urlsplit(sys.stdin.read()); print parts.username;'

    from subprocess import Popen, PIPE
    trac_admin = Popen("trac-admin . config get trac database", stdout=PIPE, shell=True)

    import urlparse
    url_parts = urlparse.urlsplit(trac_admin.stdout.read().strip())

    my_cnf = [
        "[client]",
        "host=" + url_parts.hostname,
        "user=" + url_parts.username,
        "password=" + url_parts.password,
        ""
    ]

    tmp_file.write("\n".join(my_cnf))
    tmp_file.close()

    mysqldump = Popen(["mysqldump", "--defaults-file=" + tmp_path, url_parts.path[1:]], stdout=PIPE)

    dump_path = "mysql-db-backup.sql"
    new_dump_path = "mysql-db-backup.new.sql"
    new_dump = open(new_dump_path, "w")
    new_dump.writelines(mysqldump.stdout)
    new_dump.close()

    from subprocess import call
    dump_differs = call(["diff", "-q",
                            "-I", "^INSERT INTO `session",
                            "-I", "^-- Dump completed on",
                            dump_path, new_dump_path ])

    if dump_differs:
        os.rename(new_dump_path, dump_path)
    else:
        os.unlink(new_dump_path)

    os.unlink(tmp_path)

    return not dump_differs # normal exit (0) if new dump differes, failure (1) if it doesn't
