
# global imports
from scipy.sparse import lil_matrix, hstack

# local imports
from .functions import *
from .species import *
from .density import *
from .enzymes import *
from .processes import *
            
class RbaMatrices(object):
    """
    """
    def __init__(self, data):
        ## extract information
        # extract metabolites
        self.metabolites = []
        self.external_metabolites = []
        for m in data.metabolism.species:
            if m.boundary_condition: self.external_metabolites.append(m.id)
            else: self.metabolites.append(m.id)
        # extract reactions
        reactions = data.metabolism.reactions
        self.reactions = [r.id for r in reactions]
        self.reversibility = [r.reversible for r in reactions]
        self.S = self._extract_S(self.metabolites, reactions)
        # extract functions
        self.functions = Functions(data.parameters.functions,
                                   data.parameters.aggregates)
        # extract density constraints
        compartments = [c.id for c in data.metabolism.compartments]
        self.density = Density(data.parameters.maximal_densities,
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
        """
        self.enzymes.efficiency.set_function(function_id)

    def set_medium(self, medium):
        """
        Change medium composition.
        """
        self.enzymes.efficiency.update_import(medium)
    
    def _extract_S(self, metabolites, reactions):
        S = lil_matrix((len(metabolites),len(reactions)))
        for r_index, reaction in enumerate(reactions):
            # keep ONLY internal metabolites
            for reactant in reaction.reactants:
                try:
                    m_index = self.metabolites.index(reactant.species)
                    S[m_index, r_index] = -reactant.stoichiometry
                except ValueError:
                    # external metabolite: ignore it
                    pass
            for product in reaction.products:
                try:
                    m_index = self.metabolites.index(product.species)
                    S[m_index, r_index] = product.stoichiometry
                except ValueError:
                    # external metabolite: ignore it
                    pass
        return S
