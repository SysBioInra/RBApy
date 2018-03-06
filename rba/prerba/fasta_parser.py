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
    def __init__(self, input_file):
        self.proteins = []
        self.rnas = []
        try:
            for r in SeqIO.parse(input_file, 'fasta'):
                self._create_molecule(parse_entry(r))
        except IOError:
            raise UserWarning('Please provide file ' + input_file)
        except UserWarning as e:
            raise UserWarning(input_file + ': ' + e.msg())

    def _create_molecule(self, entry):
        if entry.set_name == 'protein':
            molecule = Protein()
            molecule.cofactors = []
            self.proteins.append(molecule)
        elif entry.set_name == 'rna':
            molecule = Rna()
            self.rnas.append(molecule)
        else:
            raise UserWarning('Unknown molecule type ' + entry.set_name)
        self._initialize_molecule(molecule, entry)

    def _initialize_molecule(self, molecule, fasta_record):
        molecule.id = fasta_record.id
        molecule.stoichiometry = fasta_record.stoichiometry
        molecule.sequence = fasta_record.sequence


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
