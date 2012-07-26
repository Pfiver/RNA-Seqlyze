import os
from setuptools import setup, find_packages

import rnaseqlyze
import rnaseqlyze.worker

setup(
    name=rnaseqlyze.worker.project_name,
    version=rnaseqlyze.__version__,
    author="Patrick Pfeifer",
    author_email="patrick@patrickpfeifer.net",
    description="RNA-seq analysis | worker daemon",
    license="Mixed",
    packages=find_packages(),
    namespace_packages = ['rnaseqlyze'],
    test_suite='nose.collector',
    include_package_data=True,
    zip_safe=False,
    data_files=[('', [
        'development.ini',
        'production.ini',
        'rna-seqlyze-service',
    ])],
    setup_requires=[
        "psutil",
        "nose >= 1.1.2",
        "docopt > 0.4.1",
        "rna-seqlyze >= 0.1",
        "distribute >= 0.6.14",
    ],
    install_requires=[
        'pyramid',
        'pyramid_tm',
        'SQLAlchemy',
        'transaction',
        'waitress',
    ],
    entry_points={
        'console_scripts': [
            'rnas-worker = rnaseqlyze.worker.daemon:main',
        ],
        'paste.app_factory': [
            'main = rnaseqlyze.worker:main',
        ],
    }
)
