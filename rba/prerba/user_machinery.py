
# python 2/3 compatibility
from __future__ import division, print_function, absolute_import

# global imports
import itertools

# local imports
from rba.prerba.fasta_parser import RbaFastaParser


class UserMachinery(object):
    def __init__(self, filename, uniprot_data):
        self._filename = filename
        parser = RbaFastaParser(filename)
        self.proteins = parser.proteins
        self.rnas = parser.rnas
        self._fill_missing_protein_information(uniprot_data)

    def _fill_missing_protein_information(self, uniprot_data):
        unknown_ids = []
        for protein in self.proteins:
            try:
                uniprot = uniprot_data.create_protein_from_uniprot_id(
                    protein.id
                )
                self._copy_protein_info(protein, uniprot)
            except KeyError:
                unknown_ids.append(protein.id)
        self._print_ids_not_found_in_uniprot(unknown_ids)

    def _copy_protein_info(self, protein, uniprot_protein):
        if protein.location is None:
            protein.location = uniprot_protein.location
        if not protein.cofactors:
            protein.cofactors = uniprot_protein.cofactors
        if not protein.stoichiometry:
            protein.stoichiometry = uniprot_protein.stoichiometry
        if not protein.sequence:
            protein.sequence = uniprot_protein.sequence

    def _print_ids_not_found_in_uniprot(self, ids):
        if ids:
            print('Warning ({}): proteins {} could not be retrieved in '
                  'Uniprot.'.format(self._filename, ', '.join(ids)))

    def protein_ids(self):
        return [p.id for p in self.proteins]

    def rna_ids(self):
        return [r.id for r in self.rnas]

    def composition(self):
        return {m.id: m.stoichiometry
                for m in itertools.chain(self.proteins, self.rnas)}
