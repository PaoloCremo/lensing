# from .functions import *  # to call a F() : inf_func.F()
# from .lut import matrix   # to call matrix: inf_func.matrix

__version__ = "unknown"
try:
    from _version import __version__
except ImportError:
    # We're running in a tree that hasn't run darcsver, and didn't come with a
    # _version.py, so we don't know what our version is. This should not happen
    # very often.
    pass

from ._version import __version__
from inf_func import lut                 # to call matrix: inf_func.lut.matrix
from inf_func import functions
from inf_func import make_lut
