"""
Module defining TargetVector class.
"""

# python 2/3 compatibility
from __future__ import division, print_function

# global imports
import numpy

class TargetVector(object):
    """
    Class computing vectors where coefficients depend on growth rate.
    """

    def __init__(self, values, known_functions, def_value=None):
        """
        Constructor.

        Args:
            values: list of numeric values or function identifiers.
            known_functions: object containing function-related information.
            def_value: default value (if input value was None).
        """
        base_values = []
        self._functions = []
        for value in values:
            val = fct = None
            if value is not None:
                try:
                    val = float(value)
                except ValueError:
                    fct = known_functions.functions.get(value, None)
            else:
                val = def_value
            if val is not None:
                base_values.append(val)
                self._functions.append(None)
            elif fct is not None:
                base_values.append(0)
                self._functions.append(fct)
            else:
                raise UserWarning('Invalid value: ' + value)
        self._base_values = numpy.array(base_values, dtype='float')

    def compute(self, mu):
        """
        Compute vector for given growth rate.

        Args:
            mu: growth_rate

        Returns:
            Vector with coefficients corresponding at given growth rate.
        """
        result = self._base_values
        for i, function in enumerate(self._functions):
            if function is not None:
                result[i] = function.evaluate(mu)
        return result
