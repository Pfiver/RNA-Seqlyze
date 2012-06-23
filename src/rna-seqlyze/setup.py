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

def get_version():
    cmd = "git describe --tags"
    from subprocess import Popen, PIPE
    git_describe = Popen(cmd.split(), stdout=PIPE, stderr=PIPE)
    version, error = git_describe.communicate()
    if git_describe.returncode:
        raise Exception("couldn't determine package version: " + error)
    return version.strip()[1:]

setup(
    name='rna-seqlyze',
    version=get_version(),
    author="Patrick Pfeifer",
    author_email="patrick@patrickpfeifer.net",
    description="RNA-seq analysis | core features",
    long_description=
        "RNA-seq analysis & sequence annotation enhancement web-application",
    license="Mixed",
    packages=find_packages(),
    test_suite='nose.collector',
    entry_points={
        'console_scripts': [
            'rnas-build = rnaseqlyze.build:main',
        ],
    },
    setup_requires=[
        "nose >= 1.1.2",            # Unit testing
        "distribute >= 0.6.14",     # unsure if really required, but won't harm
    ],
    install_requires=[
        "sphinx",                   # Apidoc
        "pyflakes",                 # Style checker
        "pyramid >= 1.3.2",         # Web framework
        "MarkupSafe >= 0.15",       # Syntax highlighting in Trac
        "MySQL-python >= 1.2.2",    # Trac database backend
    ],
)
