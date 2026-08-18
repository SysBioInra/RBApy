"""
RBA package
===========

Package computing Resource Balance Analysis.
"""
import sys

if sys.version_info >= (3, 9):
    from importlib import resources
else:
    import importlib_resources as resources

from .model import RbaModel  # noqa: F401
from .prerba import *  # noqa: F401
from .core import *  # noqa: F401
from .utils import *  # noqa: F401

from . import prerba, core, utils  # noqa: F401

_version_text = resources.files(__name__).joinpath('_version.py').read_text(encoding='utf-8')
_authors_text = resources.files(__name__).joinpath('_authors.py').read_text(encoding='utf-8')

__version__ = _version_text.split("'")[1]
__author__ = _authors_text.split("'")[1]

__all__ = ['RbaModel']
__all__ += prerba.__all__
__all__ += core.__all__
__all__ += utils.__all__
