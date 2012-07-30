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
    description="RNA-seq analysis | worker daemon",
    license="Mixed",
    packages=find_packages(),
    namespace_packages = ['rnaseqlyze'],
    test_suite='nose.collector',
    include_package_data=True,
    data_files=[('', [
        'development.ini',
        'production.ini',
    ])],
    setup_requires=[
        'nose >= 1.1.2',
        'distribute >= 0.6.14',
        'setuptools_git >= 0.3',
    ],
    install_requires=[
        'docopt > 0.4.1',
        'psutil',
        'pyramid',
        'pyramid_tm',
        'waitress',
        'SQLAlchemy',
        'transaction',
        'rna-seqlyze == ' + version,
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
