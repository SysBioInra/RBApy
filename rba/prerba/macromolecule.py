"""Macromolecule information."""

# python 2/3 compatibility
from __future__ import division, print_function, absolute_import

# global imports

# local imports
from rba.prerba.default_data import DefaultData

_DEFAULT_AMINO_ACIDS = DefaultData().metabolites.aas


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
    comp = composition(sequence, 'ACGTU')
    comp['U'] += comp.pop('T')
    return comp


def aa_composition(sequence):
    return composition(sequence, _DEFAULT_AMINO_ACIDS)


class Macromolecule(object):
    def __init__(self):
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
    def __init__(self):
        super(Protein, self).__init__()
        self.cofactors = None

    def composition(self):
        comp = aa_composition(self.sequence)
        for cofactor in self.cofactors:
            comp[cofactor.chebi] = cofactor.stoichiometry
        return comp


class Rna(Macromolecule):
    def __init__(self):
        super(Rna, self).__init__()

    def composition(self):
        if isinstance(self.sequence, str):
            return ntp_composition(self.sequence)
        else:
            comp = ntp_composition(''.join(self.sequence))
            return {k: v/len(self.sequence) for k, v in comp.items()}
