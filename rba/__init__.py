"""
RBA package
===========

Package computing Resource Balance Analysis.
"""

from .model import RbaModel  # noqa: F401
from .prerba import *  # noqa: F401
from .core import *  # noqa: F401
from .utils import *  # noqa: F401

from . import prerba, core, utils  # noqa: F401

from pkg_resources import resource_string

__version__ = resource_string(__name__, '_version.py').decode("utf-8").split("'")[1]
__author__ = resource_string(__name__, '_authors.py').decode("utf-8").split("'")[1]

__all__ = ['RbaModel']
__all__ += prerba.__all__
__all__ += core.__all__
__all__ += utils.__all__
