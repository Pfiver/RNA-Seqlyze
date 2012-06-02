#!/usr/bin/env python
# encoding: utf-8

# this "distrubute" setup script depends on

    # "distribute" 0.6.26, by Tarek Ziadé
    # - source code: https://bitbucket.org/tarek/distribute
    # - documentation: http://pypi.python.org/pypi/distribute

    # to run this script, "distribute" must be installed, as follows:
    # 1) fetch the setup script
    #     $ curl -O https://bitbucket.org/tarek/distribute/raw/b69f072c0002/distribute_setup.py
    # 2) install "distribute"
    #     $ python distribute_setup.py --user

import os

from setuptools import setup
from setuptools import find_packages

setup(
    name='rna-seqlyze',
    version='0.1',
    author="Patrick Pfeifer",
    author_email="patrick@patrickpfeifer.net",
    description="RNA-seq analysis & sequence annotation enhancement web-application",
    url="http://git.pfeifer.ch/rna-seqlyze.git",
    license="Mixed",
    packages=find_packages(),
    test_suite='nose.collector',
    scripts=['scripts/' + name
        for name in os.listdir("scripts")],
    setup_requires=[
        "nose >= 1.0",
    ],
    install_requires=[
        "nose >= 1.0",
#       "Trac >= 0.13dev",
        "pyramid >= 1.3.2",
        "MarkupSafe >= 0.15",       # required by Trac
#       "MySQL-python >= 1.2.2",    # for the Trac wiki
        "distribute >= 0.6.14",     # unsure if really required, but won't harm
        "python-daemon >= 1.6",
    ],
)
