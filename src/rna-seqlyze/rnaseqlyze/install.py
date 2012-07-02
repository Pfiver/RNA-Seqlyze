"""
RNA-Seqlyze Install

This command builds and installs all software components
included with and required by the RNA-Seqlyze web application.

Usage:

    rnas-install
    rnas-install --prefix <dir>
    rnas-install -h|--help

Note:

    The command has must run from the top level RNA-Seqlyze source directory.

Options:

    --prefix <dir>  The option is passed on to the ./configure and install
                    scripts of the various prorgams that this command installs.
                    The effect is, that all produced executables will be
                    installed under the that directory.

                    The `PREFIX` variable in the "/etc/init.d/rnaseqlyze.sh"
                    worker daemon startup script _and_ the `prefix` variable
                    in the "/var/www/../rna-seqlyze.wsgi" script must both
                    be set to the directory specified here!

                    If not specified, defaults to `$HOME/.local`
"""

import os, re
from os import environ as env
from os.path import join, exists

from rnaseqlyze.build import parts, phases

def main():
    import docopt
    opts = decopt.optparse(__doc__)

    assert exists("src/rna-seqlyze/rnaseqlyze/__init__.py"), \
    "This command must be run from the top level RNA-Seqlyze source directory!"

    topdir = os.getpwd()
    prefix = opts['<dir>']
             or join(os.getenv("USER"), ".local")

    env["TOPDIR"] = topdir
    env["PREFIX"] = prefix
    env["BINDIR"] = prefix + "/bin"
    env["LIBDIR"] = prefix + "/lib"
    env["MACHTYPE"] = os.uname()[4]
    env["ARCH"] = re.sub('i.86', 'i386', env["MACHTYPE"])
    env["NCPUS_ONLN"] = str(os.sysconf("SC_NPROCESSORS_ONLN"))

    for part in parts:
        for phase in phases:
            part.execute(phase)

