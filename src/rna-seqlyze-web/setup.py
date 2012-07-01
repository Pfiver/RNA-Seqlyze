import os
from setuptools import setup, find_packages

import rnaseqlyze
import rnaseqlyze.web

setup(
    name=rnaseqlyze.web.project_name,
    version=rnaseqlyze.__version__,
    author="Patrick Pfeifer",
    author_email="patrick@patrickpfeifer.net",
    description="RNA-seq analysis | web frontend",
    license="Mixed",
    packages=find_packages(),
    namespace_packages = ['rnaseqlyze'],
    test_suite='nose.collector',
    include_package_data=True,
    zip_safe=False,
    data_files=[('', [
        'development.ini', 'production.ini'
    ])],
    setup_requires=[
        "nose >= 1.1.2",
        "rna-seqlyze >= 0.1",
        "distribute >= 0.6.14",
    ],
    install_requires=[
        'pyramid',
        'pyramid_tm',
        'pyramid_debugtoolbar',
        'SQLAlchemy',
        'transaction',
        'PasteScript',
        'zope.sqlalchemy',
    ],
    entry_points={
        'paste.app_factory': [
            'main = rnaseqlyze.web:main',
        ],
    },
)
