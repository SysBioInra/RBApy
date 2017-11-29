"""Module defining Process-related classes."""

# python 2/3 compatibility
from __future__ import division, print_function, absolute_import

# global imports
import numpy

# local imports
from rba.core.parameter_vector import ParameterVector
from rba.core import functions


class Processes(object):
    """
    Class computing process-related substructures.

    Attributes
    ----------
    ids : list of str
        Process identifiers.
    machinery : list of rba.core.species.Machinery
        Composition information of process machineries.
    capacity : rba.core.parameter_vector.ParameterVector)
        Right-hand side of process capacity constraints.
    capacity_signs : list of str
        Signs of capacity constraints ('E' for equality, 'L' for inequality).

    """

    def __init__(self, processes, species, parameters):
        """
        Constructor.

        Parameters
        ----------
        processes : rba.xml.RbaProcesses
            Structure containing process information.
        species : rba.core.species.Species
            Species information.
        parameters : rba.core.parameters.Parameters
            Parameter information

        """
        # extract ids
        self.ids = [p.id for p in processes]

        # extract machinery related information
        # extract machinery
        machinery = [p.machinery.machinery_composition for p in processes]
        self.machinery = species.create_machinery(machinery)
        # extract capacities
        # we use CPLEX conventions for signs E=equality, L=lower than.
        values = []
        self.capacity_signs = []
        for process in processes:
            if process.machinery.capacity.value is not None:
                values.append(process.machinery.capacity.value)
                self.capacity_signs.append('E')
            elif process.machinery.capacity.upper_bound is not None:
                values.append(process.machinery.capacity.value)
                self.capacity_signs.append('L')
            else:
                values.append('zero')
                self.capacity_signs.append('E')
        self.capacity = ParameterVector(values, parameters)
        # bounds
        nb_processes = len(self.ids)
        self.ub = functions.default_ub.value * numpy.ones(nb_processes)
        self.lb = numpy.zeros(nb_processes)
        self.f = numpy.zeros(nb_processes)
