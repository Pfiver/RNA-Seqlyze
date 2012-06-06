"""
The rna-seqlyze software consisty of several parts.  The majority of those parts
have been developped independent of this project and have been released under a
permissive license that allows them to be used in other (permissive licenced)
projects like this one.

This file defines a simple system and stores the commands necessary to build and
install those third-party components.
"""

from __future__ import print_function

import os, sys, re, shutil
from os import environ as env
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
        cmds = getattr(self, phase, None)
        if cmds == None: return
        print("#" * 80)
        print("# executing %s '%s' phase" % (self.name, phase))
        print("#")
        dev_null = file("/dev/null")
        logdir = "report/buildlogs"
        if not os.path.isdir(logdir): os.mkdir(logdir)
        logpath = logdir + "/%s-%s.log" % (self.name, phase)
        T = subprocess.Popen(["tee", logpath], stdin=subprocess.PIPE)
        try:
            import time
            log = lambda msg="": print(msg, file=T.stdin)
            log(time.asctime())
            log()
            log("\n".join("%s=%s" % nv
                for nv in filter(lambda i: i[0] in (
                    "PREFIX", "MACHTYPE", "NCPUS_ONLN"), env.iteritems())))
            log()
            if type(cmds) not in (list, tuple):
                cmds = cmds, # make it a 1-tuple
            for cmd in cmds:
                log("$ cd " + self.subdir)
                log("$ " + "\n  ".join(str(cmd).split("\n")))
                log()
                if type(cmd) == str:
                    ret = subprocess.call(cmd, shell=True, cwd=self.subdir,
                                stdin=dev_null, stdout=T.stdin, stderr=T.stdin)
                elif type(cmd) == MethodType:
                    def tgt():
                        sys.stdin = dev_null
                        sys.stdout = sys.stderr = T.stdin
                        os.chdir(self.subdir)
                        return cmd()
                    sp = multiprocessing.Process(target=tgt)
                    sp.start()
                    sp.join()
                    ret = sp.exitcode
                else:
                    raise Exception("Invalid '%s' phase command: %s" % (phase, repr(cmd)))
                log()
            log(time.asctime())
            if ret != 0:
                raise Exception("%s '%s' phase failed -- exit code %d" % (
                                    self.name, phase, ret))
        finally:
            T.stdin.close()
            T.wait()

# parts & phases
################

parts = []
phases = 'build', 'test', 'install'

class bcbb(Part):
    srcdir = "bcbb/nextgen"
    build = "python setup.py build"
    # save some time
    #test = "nosetests"
    install = "python setup.py install --prefix=$PREFIX"

class biopython(Part):
    build = "python setup.py build"
    # save some time
    #test = "python setup.py test"
    install = "python setup.py install --prefix=$PREFIX"

class bowtie2(Part):
    build = "make -j$NCPUS_ONLN"
    def install(self):
        """
        the bowtie2 install function
        exists because there is no 'install' target in
        the makefile, so the binaries need to be installed manually
        """
        import shutil
        for f in ("bowtie2" + x for x in ("", "-align", "-build", "-inspect")):
            shutil.copy(f, env["BINDIR"])
            os.chmod(env["BINDIR"] + "/" + f, 0775)

class samtools(Part):
    # the samtools 'build' command
    # is somewhat long so ncurses-dev is not required
    build = (
        "make -j$NCPUS_ONLN -C bcftools",
        "make -j$NCPUS_ONLN -C misc",
        """\
            make -j$NCPUS_ONLN SUBDIRS=. \\
            LIBCURSES= DFLAGS="-D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -D_USE_KNETFILE" \\
            AOBJS="bam_plcmd.o sam_view.o bam_rmdup.o bam_rmdupse.o bam_mate.o bam_stat.o bam_color.o \\
                   bamtk.o kaln.o bam2bcf.o bam2bcf_indel.o errmod.o sample.o cut_target.o phase.o bam2depth.o"
        """
    )
    install = "cp samtools $PREFIX/bin"

class cufflinks(Part):
    depends = samtools
    build = "./configure --prefix=$PREFIX --with-bam=$TOPDIR/samtools --with-eigen=$TOPDIR/eigen && make"
    install = "make install"

class kent(Part):
    build = "make -C src/lib"
    def install(self):
        """
        the kent install function
        was created because to install the few ulities required from that tree,
        the easiest way to do that is to run "make" with custom arguments for each one
        """
        for util in "bigWigToWig wigToBigWig".split(" "):
            if subprocess.call("make -C src/utils/" + util, shell=True) != 0:
                raise Exception("kent.install(): couldn't install '%s' util" % util)

class pysam(Part):
    build = "python setup.py build"
    # tests failing...
    #test = "cd tests; nosetests --exe"
    #test = "cd tests; ./pysam_test.py"
    install = "python setup.py install --prefix=$PREFIX"

class rna_seqlyze_cli(Part):
    srcdir = "rna-seqlyze-cli"
    build = "python setup.py build"
    test = "python setup.py test"
    install = "python setup.py develop --prefix=$PREFIX"

class rna_seqlyze_web(Part):
    srcdir = "rna-seqlyze-web"
    build = "python setup.py build"
    test = "python setup.py test"
    install = "python setup.py develop --prefix=$PREFIX"

class rna_seqlyze_worker(Part):
    srcdir = "rna-seqlyze-worker"
    build = "python setup.py build"
    test = "python setup.py test"
    install = "python setup.py develop --prefix=$PREFIX"

class sra_sdk(Part):
    # To get this to compile, I
    # 1) created a symlink src/sra_sdk/libxml2.so
    #    pointing to /usr/lib/libxml2.so.2 and added
    #    LDFLAGS=-L$PWD to avoid having to install libxml2-dev
    # 2) replaced the content of src/sra_sdk/libs/ext/Makefile
    #    with "all:" to skip unnesessary downloading of zlib and libbz2
    build = "LD_RUN_PATH=$LIBDIR make STATIC= STATICSYSLIBS= LDFLAGS=-L$PWD"
    def install(self):
        dir = "linux/pub/gcc/%(ARCH)s/bin/" % env
        for bin in os.listdir(dir):
            if not re.search(r'[0-9]$', bin):
                shutil.copy(dir + bin, env["BINDIR"])
                os.chmod(env["BINDIR"] + "/" + bin, 0775)
        dir = "linux/pub/gcc/%(ARCH)s/lib/" % env
        for lib in os.listdir(dir):
            if re.search(r'\.so\.[0-9]+$', lib):
                shutil.copy(dir + lib, env["LIBDIR"])

class tophat(Part):
    build = "./configure --prefix=$PREFIX --with-bam=$TOPDIR/samtools && make"
    install = "make install"

class trac(Part):
    build = "python setup.py build"
    # save some time
    #test = "python setup.py test"
    install = "python setup.py install --prefix=$PREFIX"

class trac_env(Part):
    def install(self):
        # need to discuss server
        # configuration with admin
        #destdir = "%(PREFIX)s/var/trac_env" % env
        #basedir = os.path.dirname(destdir)
        #if not os.path.isdir(basedir):
        #    os.mkdir(basedir)
        #shutil.copytree(".", destdir, symlinks=True)
        #print("Copied %s to %s\n" % (os.getcwd(), destdir))
        print("\n".join((
            "The following still needs to be done manually:",
            " 1) Set up a database",
            " 2) Restore the backup:",
            "    $ cd " + os.getcwd(),
            "    $ mysql -uUSERNAME -pPASSWORD DATABASE < mysql-db-backup.sql",
            " 4) Adjust the 'database' variable in the [trac] section in 'conf/trac.ini':",
            "    database = mysql://USERNAME:PASSWORD@localhost/DATABSE",
        )))

class transterm_hp(Part):
    build = "make"
    def install(self):
        prog = "transterm"
        data = "expterm.dat"
        shutil.copy(prog, env["BINDIR"])
        os.chmod(env["BINDIR"] + "/" + prog, 0775)
        shutil.copy(data, env["LIBDIR"])

# main routine
##############

def buildall(topdir, prefix):

    os.chdir(topdir)
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

# script enty point
###################

def main():
    import optparse
    parser = optparse.OptionParser()
    parser.add_option("-t", "--topdir",
                      dest="topdir", default=os.getcwd(),
                      help="the git source tree - default: current diretory")
    parser.add_option("-p", "--prefix",
                      dest="prefix", default=env["HOME"] + "/.local",
                      help="the installation prefix - default: $HOME/.local")
    options, args = parser.parse_args()
    buildall(options.topdir, options.prefix)
