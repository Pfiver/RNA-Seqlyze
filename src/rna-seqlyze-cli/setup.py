from pkgutil import iter_modules
from setuptools import setup, find_packages

import os, subprocess as sp
git_proc = sp.Popen("git describe --tags".split(), stdout=sp.PIPE)
git_out = git_proc.communicate()[0]
if git_proc.returncode:
    raise Exception()
version = git_out.strip()[1:]

setup(
    name=os.getcwd().split(os.sep)[-1],
    version=version,
    author="Patrick Pfeifer",
    author_email="patrick@patrickpfeifer.net",
    description="RNA-seq analysis | command-line interface",
    license="Mixed",
    packages=find_packages(),
    namespace_packages = ['rnaseqlyze'],
    test_suite='nose.collector',
    include_package_data=True,
    setup_requires=[
        'nose >= 1.1.2',
        'distribute >= 0.6.14',
        'setuptools_git >= 0.3',
    ],
    install_requires=[
        'docopt > 0.4.1',
        'SQLAlchemy',
        'rna-seqlyze == ' + version,
    ],
    entry_points={
        'console_scripts': [
            'rnas-%s = rnaseqlyze.cli.%s:main' % (n.replace("_", "-"), n) \
                for l, n, p in iter_modules(["rnaseqlyze/cli"]) if not p
        ]
    }
)
