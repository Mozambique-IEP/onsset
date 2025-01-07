"""
The onsset package contains the following modules:

 - onsset.py : main functions of the model
 - runner.py : runner is used to calibrate inputs and specify scenario runs
"""


from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version(__name__)
except PackageNotFoundError:
    # package is not installed
    __version__ = "unknown"

from .onsset import *
