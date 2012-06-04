#!/usr/bin/env python
# encoding: utf-8

import os
from setuptools import setup, find_packages

setup(
    name='rna-seqlyze-cli',
    version='0.1',
    author="Patrick Pfeifer",
    author_email="patrick@patrickpfeifer.net",
    description="RNA-seq analysis | command-line interface",
    license="Mixed",
    packages=find_packages(),
    namespace_packages = ['rnaseqlyze']
    test_suite='nose.collector',
    scripts=['scripts/' + name for name in os.listdir("scripts")],
    setup_requires=[
        "nose >= 1.0",
        "rna-seqlyze >= 0.1",
        "distribute >= 0.6.14",
    ],
)
