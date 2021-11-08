"""Module defining FastaEntry and FastaParser classes."""

# python 2/3 compatibility
from __future__ import division, print_function, absolute_import

# global imports
from collections import namedtuple
from Bio import SeqIO

# local imports
from rba.prerba.macromolecule import Protein, Rna

# class holding fasta entries
FastaEntry = namedtuple('FastaEntry',
                        'id name set_name stoichiometry sequence')


class RbaFastaParser(object):
    """Parse rba-formatted fasta file."""
    def __init__(self, input_file, uniprot_data=None):
        self._input_file = input_file
        self._uniprot = uniprot_data
        self.proteins = []
        self.rnas = []
        try:
            for r in SeqIO.parse(input_file, 'fasta'):
                self._create_molecule(parse_entry(r))
        except IOError:
            raise UserWarning('Please provide file ' + input_file)
        except UserWarning as e:
            raise UserWarning(input_file + ': ' + str(e))

    def _create_molecule(self, entry):
        if entry.set_name == 'protein':
            self.proteins.append(self._create_protein(entry))
        elif entry.set_name == 'rna':
            molecule = Rna()
            self._initialize_molecule(molecule, entry)
            self.rnas.append(molecule)
        else:
            raise UserWarning('Unknown molecule type ' + entry.set_name)

    def _create_protein(self, entry):
        result = Protein()
        self._initialize_molecule(result, entry)
        result.cofactors = []
        if self._uniprot:
            self._fill_missing_protein_information(result)
        return result

    def _fill_missing_protein_information(self, protein):
        try:
            self._copy_protein_info(
                protein,
                self._uniprot.create_protein_from_uniprot_id(protein.id)
            )
        except KeyError:
            if self._has_missing_info(protein):
                raise UserWarning(
                    '{}: protein {} has missing information and does not match'
                    ' a known UniProt id. Please fill in all information or '
                    'adapt identifier.'
                )

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
                  'UniProt.'.format(self._filename, ', '.join(ids)))

    def _initialize_molecule(self, molecule, fasta_record):
        molecule.id = fasta_record.id
        molecule.stoichiometry = fasta_record.stoichiometry
        molecule.sequence = fasta_record.sequence

    def _has_missing_info(self, molecule):
        return not (molecule.id and molecule.stoichiometry
                    and molecule.sequence)


def parse_entry(record):
    header = record.description
    try:
        [rba, id_, name, set_, sto] = header.split('|')
        sto = float(sto)
    except ValueError:
        invalid_header(header)
    if rba != 'rba':
        invalid_header(header)
    return FastaEntry(id_, name, set_, sto, str(record.seq))


def invalid_header(line):
    """Raise invalid header exception."""
    print('Invalid_header\n\t>{}'.format(line))
    print('Expected\n\t>rba|id|name|set_name|stoichiometry')
    raise UserWarning('Invalid input file')
