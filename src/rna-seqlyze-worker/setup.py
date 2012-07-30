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
    setup_requires=[
        'nose >= 1.1.2',
        'distribute >= 0.6.14',
        'setuptools_git >= 0.3',
        'docopt > 0.4.1',
        'psutil',
        'pyramid',
        'pyramid_tm',
        'waitress',
        'SQLAlchemy',
        'transaction',
        'rna-seqlyze == ' + rnaseqlyze.__version__,
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
