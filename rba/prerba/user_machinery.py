
# python 2/3 compatibility
from __future__ import division, print_function, absolute_import

# global imports
import itertools

# local imports
from rba.prerba.macromolecule import Protein
from rba.prerba.fasta_parser import parse_rba_fasta


class UserMachinery(object):
    def __init__(self, filename):
        self.proteins = []
        self.rnas = []
        for molecule in parse_rba_fasta(filename):
            if molecule.set_name == 'protein':
                self.proteins.append(self._create_protein(molecule))
            elif molecule.set_name == 'rna':
                self.rnas.append(molecule)
            else:
                print(filename + ': Unknown molecule type '
                      + molecule.set_name)

    def _create_protein(self, fasta_record):
        result = Protein()
        result.id = fasta_record.id
        result.stoichiometry = fasta_record.stoichiometry
        result.sequence = fasta_record.sequence
        return result

    def protein_ids(self):
        return [p.id for p in self.proteins]

    def rna_ids(self):
        return [r.id for r in self.rnas]

    def composition(self):
        return {m.id: m.stoichiometry
                for m in itertools.chain(self.proteins, self.rnas)}
