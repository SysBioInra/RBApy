
# python 2/3 compatibility
from __future__ import division, print_function, absolute_import

# global imports
import itertools

# local imports
from rba.prerba.macromolecule import Protein, Rna
from rba.prerba.fasta_parser import parse_rba_fasta


class UserMachinery(object):
    def __init__(self, filename):
        self.proteins = []
        self.rnas = []
        for molecule in parse_rba_fasta(filename):
            if molecule.set_name == 'protein':
                self.proteins.append(self._create_protein(molecule))
            elif molecule.set_name == 'rna':
                self.rnas.append(self._create_rna(molecule))
            else:
                print(filename + ': Unknown molecule type '
                      + molecule.set_name)

    def _create_protein(self, fasta_record):
        result = Protein()
        self._initialize_molecule(result, fasta_record)
        return result

    def _initialize_molecule(self, molecule, fasta_record):
        molecule.id = fasta_record.id
        molecule.stoichiometry = fasta_record.stoichiometry
        molecule.sequence = fasta_record.sequence
        molecule.cofactors = []

    def _create_rna(self, fasta_record):
        result = Rna()
        self._initialize_molecule(result, fasta_record)
        return result

    def protein_ids(self):
        return [p.id for p in self.proteins]

    def rna_ids(self):
        return [r.id for r in self.rnas]

    def composition(self):
        return {m.id: m.stoichiometry
                for m in itertools.chain(self.proteins, self.rnas)}
