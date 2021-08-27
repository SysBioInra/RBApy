"""Module defining Enzymes and EnzymeEfficiency classes."""

# python 2/3 compatibility
from __future__ import division, print_function, absolute_import

# global imports
import numpy

# local imports
from rba.core.functions import default_ub


class Enzymes(object):
    """
    Class computing enzyme-related substructures.

    Attributes
    ----------
    ids : list of str
        Identifiers of enzymes having a nonzero production cost.
    reaction_catalyzed : list of str
        Identifiers of reaction catalyzed by enzymes nonzero enzymes.
    machinery : rba.core.species.Machinery)
        Composition of enzymes.

    """

    def __init__(self, enzymes, species, reactions, parameters):
        """
        Constructor.

        Parameters
        ----------
        enzymes : XML node
            Structure containing enzyme information.
        species : rba.core.species.Species
            Chemical species information.
        reactions : list of str
            Reaction identifiers.
        parameters : rba.core.parameters.Parameters
            Parameter information.

        """
        # check that all reactions are found and keep only enzymes
        # that have a machinery
        reactions_left = reactions[:]
        nonzero_enzymes = []
        self.reaction_catalyzed = []
        for enzyme in enzymes.enzymes:
            reaction = enzyme.reaction
            reactions_left.remove(reaction)
            if not (enzyme.zero_cost or
                    enzyme.machinery_composition.is_empty()):
                nonzero_enzymes.append(enzyme)
                self.reaction_catalyzed.append(reaction)
        self.ids = [e.id for e in nonzero_enzymes]
        if reactions_left:
            print('Warning: did not find enzymes for following reactions: '
                  + ', '.join(reactions_left))

        # extract machinery information
        machinery = [e.machinery_composition for e in nonzero_enzymes]
        self.machinery = species.create_machinery(machinery)

        # extract efficiency information
        self._forward = [parameters[e.forward_efficiency]
                         for e in nonzero_enzymes]
        self._backward = [parameters[e.backward_efficiency]
                          for e in nonzero_enzymes]

        # bounds
        self.ub = default_ub.value * numpy.ones(len(self.ids))
        self.lb = numpy.zeros(len(self.ids))
        self.f = numpy.ones(len(self.ids))

    def efficiencies(self):
        """Compute efficiency for current parameter values."""
        forward = numpy.fromiter((fn.value for fn in self._forward),
                                 'float', len(self._forward))
        backward = numpy.fromiter((fn.value for fn in self._backward),
                                  'float', len(self._backward))
        return forward, backward
