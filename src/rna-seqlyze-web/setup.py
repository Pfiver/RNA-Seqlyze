#!/usr/bin/env python
# encoding: utf-8

import os
from setuptools import setup, find_packages

setup(
    name='rna-seqlyze-web',
    version='0.1',
    author="Patrick Pfeifer",
    author_email="patrick@patrickpfeifer.net",
    description="RNA-seq analysis | web frontend",
    license="Mixed",
    packages=find_packages(),
    namespace_packages = ['rnaseqlyze']
    test_suite='rnaseqlyze.web',
    include_package_data=True,
    zip_safe=False,
    setup_requires=[
        "nose >= 1.0",
        "rna-seqlyze >= 0.1",
        "distribute >= 0.6.14",
    ],
    install_requires=[
        'pyramid',
        'pyramid_tm',
        'pyramid_debugtoolbar',
        'SQLAlchemy',
        'transaction',
        'zope.sqlalchemy',
    ],
    entry_points={
        'paste.app_factory': [
            'main = rnaseqlyze.web:main',
        ],
        'console_scripts': [
            'initialize_rna-seqlyze-web_db = rnaseqlyze.web.scripts.initializedb:main'
        ]
    },
)
