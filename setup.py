from setuptools import setup

import os

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()


def version():
    setupDir = os.path.dirname(os.path.realpath(__file__))
    versionFile = open(os.path.join(setupDir, 'biolib', 'VERSION'))
    return versionFile.readline().strip()


setup(
    name='biolib',
    version=version(),
    author='Donovan Parks',
    author_email='donovan.parks@gmail.com',
    packages=['biolib', 'biolib.misc', 'biolib.external', 'biolib.plots'],
    package_data={'biolib': ['VERSION']},
    url='http://pypi.python.org/pypi/biolib/',
    license='GPL3',
    description='Package for common tasks in bioinformatic.',
    long_description=long_description,
    long_description_content_type="text/markdown",
    classifiers=[
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    install_requires=[],
)
