from distutils.core import setup

import os


def version():
    setupDir = os.path.dirname(os.path.realpath(__file__))
    versionFile = open(os.path.join(setupDir, 'biolib', 'VERSION'))
    return versionFile.read().strip()

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
    long_description=open('README.md').read(),
    install_requires=[
          "numpy >= 1.8.0",
          "dendropy >= 4.0.0"],
)
