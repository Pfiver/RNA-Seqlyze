#!/usr/bin/python
# encoding: utf-8

import os, sys
from os import chdir
from os.path import basename, abspath
from types import MethodType
import subprocess, multiprocessing

# a bit of infrastructure
#########################

parts = []
class PartType(type):
    def stock(cls, *ign):
        parts.append(cls())

class Part(object):

    __metaclass__ = PartType

    def __init__(self):
        self.name = self.__class__.__name__
        try:
            self.subdir = "src/" + self.srcdir
        except AttributeError:
            self.subdir = "src/" + self.name

    def run(self, command):
        try:
            cmd = getattr(self, command)
        except AttributeError:
            return
        print "running %s '%s' command" % (self.name, command)
        dev_null = file("/dev/null")
        tee_log = subprocess.Popen(
                    ["tee", "report/buildlogs/%s-%s.log" % (
                                self.name, command)], stdin=subprocess.PIPE)
        try:
            if type(cmd) == str:
                ret = subprocess.call(cmd, shell=True, cwd=self.subdir,
                                      stdin=dev_null,
                                      stdout=tee_log.stdin, stderr=tee_log.stdin)
            elif type(cmd) == MethodType:
                def target():
                    sys.stdin = dev_null
                    sys.stdout = sys.stderr = tee_log.stdin
                    os.chdir(self.subdir)
                    return cmd()
                sp = multiprocessing.Process(target=target)
                sp.start()
                sp.join()
                ret = sp.exitcode
            else:
                raise Exception("Invalid %s command: %s" % (command, repr(cmd)))
            if ret != 0:
                raise Exception("%s '%s' command failed -- exit code %d" % (
                                    self.name, command, ret))
        finally:
            tee_log.stdin.close()
            tee_log.wait()

# start auto-creating and stocking
# "Part" instances upon creation of (sub)classes
PartType.__init__ = PartType.stock

# the parts
###########

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

    chdir(topdir)
    os.environ["PREFIX"] = prefix
    os.environ["MACHTYPE"] = os.uname()[4]
    os.environ["NCPUS_ONLN"] = str(os.sysconf("SC_NPROCESSORS_ONLN"))

    for command in 'build', 'test', 'install':
        for part in parts:
            part.run(command)
