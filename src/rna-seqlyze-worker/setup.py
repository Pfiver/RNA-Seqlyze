import os
from setuptools import setup, find_packages

setup(
    name='rna-seqlyze-worker',
    version='0.1',
    author="Patrick Pfeifer",
    author_email="patrick@patrickpfeifer.net",
    description="RNA-seq analysis | worker daemon",
    license="Mixed",
    packages=find_packages(),
    namespace_packages = ['rnaseqlyze'],
    test_suite='rnaseqlyze.worker',
    include_package_data=True,
    zip_safe=False,
    setup_requires=[
        "nose >= 1.1.2",
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
        'paste.app_factory': [
            'main = rnaseqlyze.worker:wsgi_app',
        ],
    }
)
