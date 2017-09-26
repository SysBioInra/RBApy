"""Module defining RbaMatrices class."""

# python 2/3 compatibility
from __future__ import division, print_function

# local imports
from rba.core.metabolism import Metabolism
from rba.core.functions import Functions
from rba.core.species import Species
from rba.core.density import Density
from rba.core.enzymes import Enzymes
from rba.core.processes import Processes


class RbaMatrices(object):
    """
    Class converting RBA data into computing substructures.

    Attributes:
        metabolites: list of internal metabolites.
        external_metabolites: list of external metabolites.
        reactions: list of reaction identifiers.
        reversibility: list indicating whether reactions are reversible.
        S: stoichiometry matrix
        functions (functions.Functions): substructure storing functions.
        density (density.Density): substructure storing density information.
        species (species.Species): substructure storing species-related
            information.
        enzymes (enzymes.Enzymes): substructure storing enzyme-related
            information.
        processes (processes.Processes): substructure storing processes.

    """

    def __init__(self, data):
        """
        Constructor.

        Args:
            data (rba.RbaModel): RBA model containing raw data.
        """
        # extract metabolism
        self.metabolism = Metabolism(data.metabolism.species,
                                     data.metabolism.reactions)
        # extract functions
        self.functions = Functions(data.parameters.functions,
                                   data.parameters.aggregates)
        # extract density constraints
        compartments = [c.id for c in data.metabolism.compartments]
        self.density = Density(data.parameters.target_densities,
                               self.functions, compartments)
        # extract base species composition (metabolites + polymers)
        self.species = Species(data, self.metabolism.internal)
        # extract enzyme information
        self.enzymes = Enzymes(data.enzymes, self.species,
                               self.metabolism.reactions)
        # add synthesis reaction for metabolites that are also macromolecules
        (new_reactions, names) = self.species.metabolite_synthesis()
        nb_reactions = len(new_reactions)
        if nb_reactions > 0:
            self.metabolism.add_reactions(new_reactions, names,
                                          [False] * nb_reactions)
        # extract process information
        self.processes = Processes(data.processes.processes,
                                   self.species, self.functions)
        # setup medium
        self.set_medium(data.medium)

    def set_catalytic_function(self, function_id):
        """
        Change set of catalytic functions to use.

        Args:
            function_id: identifier matching an efficiency function defined
                in the enzyme substructure.
        """
        self.enzymes.efficiency.set_function(function_id)

    def set_medium(self, medium):
        """
        Change medium composition.

        Args:
            medium: dict mapping metabolite prefixes with their concentration.
        """
        self.metabolism.set_boundary_fluxes(medium)
        self.enzymes.efficiency.update_import(medium)
