"""Module defining Density class."""

# python 2/3 compatibility
from __future__ import division, print_function, absolute_import

# local imports
from rba.core.parameter_vector import ParameterVector


class Density(object):
    """
    Class computing density-related substructures.

    Attributes
    ----------
    compartments : list of str
        compartment identifiers. compartments[i] is
        the identifier of the compartment involved in the ith constraint.
    compartment_indices : list of int
        compartment_indices[i] is the index (in the list of ids provided at
        construction) of the compartment involved in the ith constraint.
    signs : list of str
        signs of constraints ('E' for equality, 'L' for
        inequality - Lower than)
    values : rba.core.parameter_vector.ParameterVector
        Right-hand side of constraints (depending on μ).

    """

    def __init__(self, target_densities, parameters, known_compartments):
        """
        Constructor.

        Parameters
        ----------
        target_densities : rba.xml.ListOfTargetDensities
            Structure holding target density information.
        parameters : rba.core.parameters.Parameters
            Parameter information.
        known_compartments : list of str
            ids of compartments in the system.

        """
        # extract target densities
        self.compartments = [md.compartment for md in target_densities]
        self.compartment_indices = [known_compartments.index(md.compartment)
                                    for md in target_densities]
        values = []
        self.signs = []
        for target in target_densities:
            if target.value is not None:
                values.append(target.value)
                self.signs.append('E')
            elif target.upper_bound is not None:
                values.append(target.upper_bound)
                self.signs.append('L')
            else:
                raise UserWarning('Density constraint ' + target.compartment
                                  + ': you must specify a value or an upper '
                                  'bound.')
        self.values = ParameterVector(values, parameters)
