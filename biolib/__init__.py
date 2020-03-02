import os


def version():
    """Read program version from file."""
    binDir = os.path.dirname(os.path.realpath(__file__))
    versionFile = open(os.path.join(binDir, 'VERSION'))
    return versionFile.readline().strip()

__version__ = version()
