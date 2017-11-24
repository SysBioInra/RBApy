"""Module defining Metabolism class."""

# python 2/3 compatibility
from __future__ import division, print_function, absolute_import

# global imports
from scipy.sparse import lil_matrix, hstack
from itertools import chain


class Metabolism(object):
    """
    Class computing metabolism-related substructures.

    Attributes
    ----------
    internal : list of str
        Identifiers of internal metabolites.
    external : list of str
        Identifiers of external metabolites.
    reactions: list of str
        Identifiers of reactions.
    reversibility: list of bool
        Reversibility of reactions.
    S : sparse matrix
        Stoichiometry matrix (with only internal metabolites).
    zero_lb : list
        Indices of reactions with lower bound forced to zero because of
        medium.
    zero_ub : list
        Indices of reactions with upper bound forced to zero because of
        medium.

    """

    def __init__(self, metabolites, reactions):
        """
        Constructor.

        Parameters
        ----------
        metabolites: rba.xml.ListOfSpecies
            Data structure with metabolite information.
        reactions: rba.xml.ListOfReactions
            Data structure with metabolite information.

        """
        self.internal = [m.id for m in metabolites if not m.boundary_condition]
        self.external = [m.id for m in metabolites if m.boundary_condition]
        self.reactions = [r.id for r in reactions]
        self.reversibility = [r.reversible for r in reactions]
        full_S = build_S(self.external + self.internal, reactions)
        S_ext = full_S[:len(self.external), ]
        self.S = full_S[len(self.external):, ]
        # find boundary reactions (involving external metabolites)
        self._boundary_lb = {}
        self._boundary_ub = {}
        for met_id, r_indices, values \
                in zip(self.external, S_ext.rows, S_ext.data):
            met_prefix = met_id.rsplit('_', 1)[0]
            self._boundary_lb[met_prefix] = [
                i for i, v in zip(r_indices, values) if v > 0
                ]
            self._boundary_ub[met_prefix] = [
                i for i, v in zip(r_indices, values) if v < 0
                ]
        self.zero_lb = []
        self.zero_ub = []

    def set_boundary_fluxes(self, medium):
        """
        Find zero lower and upper bounds imposed by external medium.

        Parameters
        ----------
        medium : dict
            Dictionary where keys are metabolite prefixes and values are
            concentrations in external medium.

        """
        zero_metabolites = [m for m, c in medium.items() if c == 0]
        self.zero_lb = list(set(
            chain.from_iterable(self._boundary_lb[m] for m in zero_metabolites)
            ))
        self.zero_ub = list(set(
            chain.from_iterable(self._boundary_ub[m] for m in zero_metabolites)
            ))

    def add_reactions(self, reactions, names, reversibility):
        """
        Add reactions.

        Parameters
        ----------
        reactions : list of vectors
            Reactions in vector format where order must match order of
            internal metabolites.
        names : list of str
            Identifiers of reactions.
        reversibility : list of bool
            Reversibility of reactions.

        """
        assert(len(reactions) == len(names) == len(reversibility))
        self.S = hstack([self.S] + reactions)
        self.reactions += names
        self.reversibility += reversibility


def build_S(metabolites, reactions):
    """
    Build stoichiometry matrix from metabolites and reactions.

    Parameters
    ----------
    metabolites:
        Metabolite identifiers (used to define row order).
    reactions: rba.xml.ListOfReactions
        Reaction data.

    Returns
    -------
    scipy.sparse.lil_matrix
        Stoichiometry matrix.

    """
    m_index = {m: i for i, m in enumerate(metabolites)}
    S = lil_matrix((len(metabolites), len(reactions)))
    for r_index, reaction in enumerate(reactions):
        for reactant in reaction.reactants:
            S[m_index[reactant.species], r_index] = -reactant.stoichiometry
        for product in reaction.products:
            S[m_index[product.species], r_index] = product.stoichiometry
    return S
