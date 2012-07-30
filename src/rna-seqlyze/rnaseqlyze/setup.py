"""
RNA-Seqlyze Setup

This command builds and installs all 3rd-party software
components included with and required by the RNA-Seqlyze application.

Usage:
    rnas-setup
    rnas-setup --prefix <dir>
    rnas-setup -h|--help

Note:
    The command has must run from the top level RNA-Seqlyze source directory.

Options:
    --prefix <dir>
                    The option is passed on to the ./configure and install
                    scripts of the various prorgams that this command installs.
                    The effect is, that all produced executables will be
                    installed under the that directory.

                    If not specified, defaults to ``$HOME/.local``
"""

import os, re
from os import environ as env
from os.path import join, exists

from rnaseqlyze.build import parts, phases

def main():
    import docopt
    opts = docopt.docopt(__doc__)

    assert exists("src/rna-seqlyze/rnaseqlyze/__init__.py"), \
    "This command must be run from the top level RNA-Seqlyze source directory!"

    topdir = os.getcwd()
    prefix = opts['--prefix'] \
             or join(os.getenv("HOME"), ".local")

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

    print """\

 All 3rd-party software has been sucessfully installed under

    PREFIX=%s

""" % prefix

