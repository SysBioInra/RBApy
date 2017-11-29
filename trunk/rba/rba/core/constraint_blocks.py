"""Module defining ConstraintBlocks class."""

# python 2/3 compatibility
from __future__ import division, print_function, absolute_import

# local imports
from rba.core.metabolism import Metabolism
from rba.core.parameters import Parameters
from rba.core.species import Species
from rba.core.density import Density
from rba.core.enzymes import Enzymes
from rba.core.processes import Processes
from rba.core.targets import Targets


class ConstraintBlocks(object):
    """
    Class converting RBA data into computing substructures.

    Attributes
    ----------
    metabolism : rba.core.metabolism.Metabolism
        Metabolism information.
    parameters : rba.core.parametres.Parameters
        Parameters.
    density : rba.core.density.Density
        Density information.
    species : rba.core.species.Species)
        Species information.
    enzymes : rba.core.enzymes.Enzymes
        Enzyme information.
    processes : rba.core.processes.Processes
        Process information.
    targets : rba.core.targets.Targets
        Target information.

    """

    def __init__(self, data):
        """
        Constructor.

        Parameters
        ----------
        data : rba.RbaModel
            RBA model containing raw data.

        """
        # extract metabolism
        self.metabolism = Metabolism(data.metabolism.species,
                                     data.metabolism.reactions)
        # extract parameters
        self.parameters = Parameters(data.parameters.functions,
                                     data.parameters.aggregates)
        # extract density constraints
        compartments = [c.id for c in data.metabolism.compartments]
        self.density = Density(data.density.target_densities,
                               self.parameters, compartments)
        # extract base species composition (metabolites + polymers)
        self.species = Species(data, self.metabolism.internal)
        # extract enzyme information
        self.enzymes = Enzymes(data.enzymes, self.species,
                               self.metabolism.reactions, self.parameters)
        # add synthesis reaction for metabolites that are also macromolecules
        (new_reactions, names) = self.species.metabolite_synthesis()
        nb_reactions = len(new_reactions)
        if nb_reactions > 0:
            self.metabolism.add_reactions(new_reactions, names,
                                          [False] * nb_reactions)
        # extract process information
        self.processes = Processes(data.processes.processes,
                                   self.species, self.parameters)
        # extract target information
        self.targets = Targets(data.targets,
                               self.species, self.parameters)
        # setup medium
        self.set_medium(data.medium)

    def set_medium(self, medium):
        """
        Change medium composition.

        Args:
            medium: dict mapping metabolite prefixes with their concentration.
        """
        self.metabolism.set_boundary_fluxes(medium)
        self.parameters.update_medium(medium)
