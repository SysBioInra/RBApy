
# python 2/3 compatibility
from __future__ import division, print_function, absolute_import

# global imports
import itertools

# local imports
from rba.prerba.fasta_parser import RbaFastaParser


class UserMachinery(object):
    def __init__(self, filename, uniprot_data):
        parser = RbaFastaParser(filename, uniprot_data)
        self.proteins = parser.proteins
        self.rnas = parser.rnas

    def protein_ids(self):
        return [p.id for p in self.proteins]

    def rna_ids(self):
        return [r.id for r in self.rnas]

    def composition(self):
        return {m.id: m.stoichiometry
                for m in itertools.chain(self.proteins, self.rnas)}
