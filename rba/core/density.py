"""
Module defining Density class.
"""

# python 2/3 compatibility
from __future__ import division, print_function

# local imports
from rba.core.target_vector import TargetVector

class Density(object):
    """
    Class computing density-related substructures.

    Attributes:
        compartments: list of compartment identifiers. compartments[i] is
            the identifier of the compartment involved in the ith constraint.
        compartment_indices: List of indices. compartment_indices[i] is
            the index (in the list of ids provided at construction) of
            the compartment involved in the ith constraint.
        signs: list of signs of constraints ('E' for equality, 'L' for
            inequality - Lower than)
        values (target_vector.TargetVector): object used to compute the
            right-hand side of constraints depending on mu.
    """

    def __init__(self, target_densities, known_functions, known_compartments):
        """
        Constructor.

        Args:
            target_densities: xml structure holding target density information.
            known_functions: dict mapping function ids with the object used to
                compute them.
            known_compartments: list of ids of compartments in the system.
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
        self.values = TargetVector(values, known_functions)
