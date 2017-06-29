"""
Module defining RbaMatrices class.
"""

# python 2/3 compatibility
from __future__ import division, print_function

# global imports
from scipy.sparse import lil_matrix, hstack

# local imports
from rba.rba_main.functions import Functions
from rba.rba_main.species import Species
from rba.rba_main.density import Density
from rba.rba_main.enzymes import Enzymes
from rba.rba_main.processes import Processes

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
        ## extract information
        # extract metabolites
        self.metabolites = []
        self.external_metabolites = []
        for species in data.metabolism.species:
            if species.boundary_condition:
                self.external_metabolites.append(species.id)
            else:
                self.metabolites.append(species.id)
        # extract reactions
        reactions = data.metabolism.reactions
        self.reactions = [r.id for r in reactions]
        self.reversibility = [r.reversible for r in reactions]
        self.S = build_stoichiometry_matrix(self.metabolites, reactions)
        # extract functions
        self.functions = Functions(data.parameters.functions,
                                   data.parameters.aggregates)
        # extract density constraints
        compartments = [c.id for c in data.metabolism.compartments]
        self.density = Density(data.parameters.target_densities,
                               self.functions, compartments)
        # extract base species composition (metabolites + polymers)
        self.species = Species(data, self.metabolites)
        # extract enzyme information
        self.enzymes = Enzymes(data.enzymes, self.species, self.reactions)
        # add synthesis reaction for metabolites that are also macromolecules
        (new_reactions, names) = self.species.metabolite_synthesis()
        if len(new_reactions) > 0:
            nb_reactions = len(new_reactions)
            self.S = hstack([self.S] + new_reactions)
            self.reactions += names
            self.reversibility += [False] * nb_reactions
        # extract process information
        self.processes = Processes(data.processes.processes,
                                   self.species, self.functions)

        ## setup medium
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
        self.enzymes.efficiency.update_import(medium)

def build_stoichiometry_matrix(metabolites, reactions):
    """
    Build stoichiometry matrix from metabolites and reactions.

    Args:
        metabolites: list of internal metabolites.
        reactions: RBA reaction data.
    
    Returns:
        Stoichiometry matrix in sparse format (excluding external metabolites).
    """
    S = lil_matrix((len(metabolites), len(reactions)))
    for r_index, reaction in enumerate(reactions):
        # keep ONLY internal metabolites
        for reactant in reaction.reactants:
            try:
                m_index = metabolites.index(reactant.species)
                S[m_index, r_index] = -reactant.stoichiometry
            except ValueError:
                # external metabolite: ignore it
                pass
        for product in reaction.products:
            try:
                m_index = metabolites.index(product.species)
                S[m_index, r_index] = product.stoichiometry
            except ValueError:
                # external metabolite: ignore it
                pass
    return S
