#!/usr/bin/python
# encoding: utf-8

"""
The rna-seqlyze software consisty of several parts.  The majority of those parts
have been developped independent of this project and have been released under a
permissive license that allows them to be used in other (permissive licenced)
projects like this one.

This file defines a simple system and stores the commands necessary to build and
install those third-party components.
"""

from __future__ import print_function

import os, sys
from types import MethodType
import subprocess, multiprocessing

# a bit of infrastructure
#########################

class PartType(type):

    def __init__(cls, *ign):
        """
        auto-create and stock instances
        upon creation of "Part" (sub)classes

        appends a new instance of the
        created class to the "parts" list, if it exists
        """
        try:
            parts.append(cls())
        except NameError:
            pass

class Part(object):

    __metaclass__ = PartType

    def __init__(self):
        self.name = self.__class__.__name__
        try:
            self.subdir = "src/" + self.srcdir
        except AttributeError:
            self.subdir = "src/" + self.name

    def execute(self, phase):
        try:
            cmd = getattr(self, phase)
        except AttributeError:
            return
        print("#" * 80)
        print("# executing %s '%s' phase" % (self.name, phase))
        print("#")
        dev_null = file("/dev/null")
        tee_log = subprocess.Popen(
                    ["tee", "report/buildlogs/%s-%s.log" % (
                                self.name, phase)], stdin=subprocess.PIPE)
        try:
            import time
            envars = "PREFIX", "MACHTYPE", "NCPUS_ONLN"
            env = filter(lambda i: i[0] in envars, os.environ.iteritems())
            log = lambda msg="": print(msg, file=tee_log.stdin)
            log(time.asctime())
            log()
            log("\n".join("%s=%s" % nv for nv in env))
            log()
            log("$ " + "\n  ".join(str(cmd).split("\n")))
            log()
            if type(cmd) == str:
                ret = subprocess.call(cmd, shell=True,
                        cwd=self.subdir, stdin=dev_null,
                        stdout=tee_log.stdin, stderr=tee_log.stdin)
            elif type(cmd) == MethodType:
                def tgt():
                    sys.stdin = dev_null
                    sys.stdout = sys.stderr = tee_log.stdin
                    os.chdir(self.subdir)
                    return cmd()
                sp = multiprocessing.Process(target=tgt)
                sp.start()
                sp.join()
                ret = sp.exitcode
            else:
                raise Exception("Invalid %s phase: %s" % (phase, repr(cmd)))
            log()
            log(time.asctime())
            if ret != 0:
                raise Exception("%s '%s' phase failed -- exit code %d" % (
                                    self.name, phase, ret))
        finally:
            tee_log.stdin.close()
            tee_log.wait()

# parts & phases
################

parts = []
phases = 'build', 'test', 'install'

class bcbb(Part):
    srcdir = "bcbb/nextgen"
    build = "python setup.py build"
#    test = "nosetests"
#    install = "python setup.py install --user"

class biopython(Part):
    build = "python setup.py build"
#    test = "python setup.py test"
    install = "python setup.py install --user"

class bowtie2(Part):
    build = "make -j$NCPUS_ONLN"
    def install(self):
        """
        the bowtie2 install function
        exists because there is no 'install' target in
        the makefile, so the binaries need to be installed manually
        """
        import shutil
        dst = os.environ["PREFIX"] + "/bin"
        for f in ("bowtie2" + x for x in ("", "-align", "-build", "-inspect")):
            shutil.copy(f, dst)
            os.chmod(os.path.join(dst, f), 0775)

class samtools(Part):
    # the samtools 'build' command
    # is somewhat long so ncurses-dev is not required
    build = """\
        make -j$NCPUS_ONLN -C bcftools &&
        make -j$NCPUS_ONLN -C misc &&
        make -j$NCPUS_ONLN SUBDIRS=. \
            LIBCURSES= DFLAGS="-D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -D_USE_KNETFILE" \
            AOBJS="bam_plcmd.o sam_view.o bam_rmdup.o bam_rmdupse.o bam_mate.o bam_stat.o bam_color.o \
                   bamtk.o kaln.o bam2bcf.o bam2bcf_indel.o errmod.o sample.o cut_target.o phase.o bam2depth.o"
    """
    install = "cp samtools $PREFIX/bin"

class cufflinks(Part):
    depends = samtools
    build = "./configure --prefix=$PREFIX --with-bam=$PWD/../samtools --with-eigen=$PWD/../eigen && make"
    install = "make install"

class kent(Part):
    build = "make -C src/lib"
    def install(self):
        """
        the kent utils install function
        was created because to install those it is
        easiest to employ make with some boilerplate arguments
        """
        for util in "bigWigToWig wigToBigWig".split(" "):
            if subprocess.call("make BINDIR=$PREFIX/bin -C src/utils/" + util, shell=True) != 0:
                raise Exception("kent.install(): couldn't install '%s' util" % util)

# main routine
##############

def buildall(topdir, prefix):

    from os import chdir

    chdir(topdir)
    os.environ["PREFIX"] = prefix
    os.environ["MACHTYPE"] = os.uname()[4]
    os.environ["NCPUS_ONLN"] = str(os.sysconf("SC_NPROCESSORS_ONLN"))

    for part in parts:
        for phase in phases:
            part.execute(phase)
