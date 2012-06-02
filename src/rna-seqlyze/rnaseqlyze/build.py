#!/usr/bin/python
# encoding: utf-8

import os
from os import chdir
from os.path import basename, abspath
from types import FunctionType
import subprocess, multiprocessing

class Part(object):

    def __init__(self, srcdir, depends=None, **commands):
        self.name = srcdir.replace("/", "_")
        self.subdir = "src/" + srcdir
        self.depends = depends
        self.commands = commands

    def run(self, command):
        if command not in self.commands:
            return
        devnull = file("/dev/null")
        logfile = file("report/buildlogs/%s-%s.log" % (self.name, command), "w")
        cmd = self.commands[command]
        if type(cmd) == 'str':
            subprocess.call(cmd, shell=True, cwd=subdir,
                            stdin=devnull, stdout=logfile, stderr=logfile)
        elif type(cmd) == FunctionType:
            sp = multiprocessing.Process(target=cmd, args=(self,))
            sp.start()
            sp.join()

def bowtie2_install(part):
    print("bowtie2_install(part): part.subdir = '%s'" % part.subdir)

parts = [

    Part(srcdir="bcbio/nextgen",
        build="python setup.py build",
        test="nosetests",
        install="python setup.py install --user"),

    Part(srcdir="biopython",
        build="python setup.py build",
        test="python setup.py test",
        install="python setup.py install"),

    Part(srcdir="bowtie2",
        build="make",
        install=bowtie2_install),

    Part(srcdir="samtools",
        build="make",
        install="cp samtools $PREFIX/bin"),

    Part(srcdir="cufflinks",
        depends=["samtools"],
        build="./configure --prefix=$PREFIX && make",
        install="make install"),
]

def buildall(topdir, prefix):

    chdir(topdir)
    os.environ["PREFIX"] = prefix

    for command in 'build', 'test', 'install':
        for part in parts:
            part.run(command)
