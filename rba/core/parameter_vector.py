"""Module defining TargetVector class."""

# python 2/3 compatibility
from __future__ import division, print_function, absolute_import

# global imports
import numpy

# local imports
# from rba.core.functions


class ParameterVector(object):
    """Vectors with coefficients depending on growth rate."""

    def __init__(self, parameter_list, known_parameters):
        """
        Constructor.

        Parameters
        ----------
        parameter_list : list of str
            Function identifiers used to compute the vector.
        known_parameters : rba.core.functions.Parameters
            Parameter information.
        """
        self._functions = []
        for parameter in parameter_list:
            self._functions.append(known_parameters[parameter])

    def compute(self):
        """Compute vector with current parameter values."""
        return numpy.fromiter((fn.value for fn in self._functions),
                              'float', len(self._functions))
