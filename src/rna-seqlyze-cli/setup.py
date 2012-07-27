from pkgutil import iter_modules
from setuptools import setup, find_packages

import rnaseqlyze
import rnaseqlyze.cli

setup(
    name=rnaseqlyze.cli.project_name,
    version=rnaseqlyze.__version__,
    author="Patrick Pfeifer",
    author_email="patrick@patrickpfeifer.net",
    description="RNA-seq analysis | command-line interface",
    license="Mixed",
    packages=find_packages(),
    namespace_packages = ['rnaseqlyze'],
    test_suite='nose.collector',
    setup_requires=[
        "SQLAlchemy",
        "nose >= 1.1.2",
        "docopt > 0.4.1",
        "pyinotify >= 0.9",
        "rna-seqlyze >= 0.1",
        "distribute >= 0.6.14",
    ],
    entry_points={
        'console_scripts': [
            'rnas-%s = rnaseqlyze.cli.%s:main' % (n.replace("_", "-"), n) \
                for l, n, p in iter_modules(["rnaseqlyze/cli"]) if not p
        ]
    }
)
