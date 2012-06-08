# this "distrubute" setup script depends on

    # "distribute" 0.6.26
    # - source code: https://bitbucket.org/tarek/distribute
    # - documentation: http://pypi.python.org/pypi/distribute

    # to run this script, "distribute" must be installed, as follows:
    # 1) fetch the setup script
    #     $ curl -O https://bitbucket.org/tarek/distribute/raw/b69f072c0002/distribute_setup.py
    # 2) install "distribute"
    #     $ python distribute_setup.py --user

import os
from setuptools import setup, find_packages

setup(
    name='rna-seqlyze',
    version='0.1',
    author="Patrick Pfeifer",
    author_email="patrick@patrickpfeifer.net",
    description="RNA-seq analysis | core features",
    long_description="RNA-seq analysis & sequence annotation enhancement web-application",
    license="Mixed",
    packages=find_packages(),
    test_suite='nose.collector',
    entry_points={
        'console_scripts': [
            'build-rnaseqlyze = rnaseqlyze.build:main',
            'rnaseqlyze-dbinit = rnaseqlyze.dbinit:main',
        ],
    },
    setup_requires=[
        "nose >= 1.0",
        "distribute >= 0.6.14",     # unsure if really required, but won't harm
    ],
    install_requires=[
        "pyflakes",                 # Buildbot
        "nose >= 1.1.2",
        "pyramid >= 1.3.2",
        "MarkupSafe >= 0.15",       # in Trac
        "MySQL-python >= 1.2.2",    # in Trac
    ],
)
