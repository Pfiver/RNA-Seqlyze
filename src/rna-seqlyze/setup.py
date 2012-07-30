# rna-seqlyze setup script

    # depends on git and "distribute"
    # - git: git-scm.com
    # - "distribute":
    #  - docs: http://packages.python.org/distribute/
    #  - code: https://bitbucket.org/tarek/distribute
    #  - installation:
    #    $ curl -O http://python-distribute.org/distribute_setup.py
    #    $ python distribute_setup.py --user

from setuptools import setup, find_packages

import os, subprocess as sp
git_proc = sp.Popen("git describe --tags".split(), stdout=sp.PIPE)
git_out = git_proc.communicate()[0]
if git_proc.returncode:
    raise Exception()

setup(
    name=os.getcwd().split(os.sep)[-1],
    version=git_out.strip()[1:],
    author="Patrick Pfeifer",
    author_email="patrick@patrickpfeifer.net",
    url="http://biocalc.fhnw.ch/",
    platforms=["linux2"],
    description="RNA-seq analysis | core features",
    long_description=
        "RNA-seq analysis & sequence annotation enhancement web-application",
    license="Mixed",
    packages=find_packages(),
    test_suite='nose.collector',
    include_package_data=True,
    data_files=[('', [
        'rnaseqlyze.ini',
    ])],
    setup_requires=[
        'nose >= 1.1.2',            # Unit testing
        'distribute >= 0.6.14',     # Packaging (setuptools)
        'setuptools_git >= 0.3',    # Packaging (include_package_data)
    ],
    install_requires=[
        'sphinx',                   # Apidoc
        'pyflakes',                 # Style checker
        'SQLAlchemy',
        'docopt > 0.4.1',
        'pyramid >= 1.3.2',         # Web framework
        'MarkupSafe >= 0.15',       # Syntax highlighting in Trac
    ],
    entry_points={
        'console_scripts': [
            'rnas-setup = rnaseqlyze.setup:main',
        ],
    },
)
