"""Macromolecule information."""

# python 2/3 compatibility
from __future__ import division, print_function, absolute_import

# global imports

# local imports
from rba.prerba.default_data import DefaultData


class Macromolecule(object):
    def __init__(self):
        """Build default protein."""
        self.id = None
        self.location = None
        self.stoichiometry = None
        self.sequence = None


class Protein(Macromolecule):
    """
    Basic protein information.

    Attributes
    ----------
    location : str
        Location of protein.
    stoichiometry : int or float
        Stoichiometry of protein within enzymatic complex.
    cofactors : Cofactor
        Cofactors of protein.
    sequence : str
        Sequence of protein (amino acids in one-lette format).

    """
    DEFAULT_AMINO_ACIDS = DefaultData().metabolites.aas

    def __init__(self):
        """Build default protein."""
        self.cofactors = []

    def composition(self):
        comp = self._aa_composition()
        for cofactor in self.cofactors:
            comp[cofactor.chebi] = cofactor.stoichiometry
        return comp

    def _aa_composition(self):
        return composition(self.sequence, self.DEFAULT_AMINO_ACIDS)


def composition(sequence, alphabet):
    """Compute composition of sequence with given alphabet."""
    comp = dict.fromkeys(alphabet, 0)
    for n in sequence:
        try:
            comp[n] += 1
        except KeyError:
            pass
    return comp


def ntp_composition(sequence):
    """Translate rna or dna sequence into ntp composition."""
    comp = composition(sequence, 'ACGTU')
    comp['U'] += comp.pop('T')
    return comp
