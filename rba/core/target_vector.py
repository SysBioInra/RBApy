"""Module defining TargetVector class."""

# python 2/3 compatibility
from __future__ import division, print_function, absolute_import

# global imports
import numpy

# local imports
from rba.core import functions


class TargetVector(object):
    """Vectors with coefficients depending on growth rate."""

    def __init__(self, parameter_list, known_parameters, def_value=None):
        """
        Constructor.

        Parameters
        ----------
        parameter_list : list of str
            Function identifiers used to compute the vector.
        known_parameters : rba.core.functions.Parameters
            Parameter information.
        def_value : float
            Default value (if some input paramater was set to None).
        """
        base_values = []
        self._functions = []
        for parameter in parameter_list:
            if parameter:
                base_values.append(0)
                self._functions.append(known_parameters[parameter])
            elif def_value is not None:
                base_values.append(def_value)
                self._functions.append(None)
            else:
                raise UserWarning('Default value is missing...')
        self._base_values = numpy.array(base_values, dtype='float')

    def compute(self):
        """Compute vector with current parameter values."""
        result = self._base_values
        for i, function in enumerate(self._functions):
            if function is not None:
                result[i] = function.value
        return result
