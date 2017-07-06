"""
Module defining FastaEntry and FastaParser classes.
"""

# python 2/3 compatibility
from __future__ import division, print_function

# global imports
from collections import namedtuple
from Bio import SeqIO

# class holding fasta entries
FastaEntry = namedtuple('FastaEntry', 'id name set_name stoichiometry sequence')

class FastaParser(object):
    """
    Module parsing rba-formatted fasta files.

    Attributes:
        entries: list of FastaEntry read.
    """

    def __init__(self, input_file):
        self.entries = []
        try:
            for record in SeqIO.parse(input_file, 'fasta'):
                header = record.description
                try:
                    [rba, id_, name, set_, sto] = header.split('|')
                    sto = float(sto)
                except ValueError:
                    self._invalid_header(header)
                if rba != 'rba':
                    self._invalid_header(header)
                self.entries.append(FastaEntry(id_, name, set_,
                                               sto, str(record.seq)))
        except IOError:
            raise UserWarning('Please provide file ' + input_file)

    @staticmethod
    def _invalid_header(line):
        print('Invalid_header\n\t>{}'.format(line))
        print('Expected\n\t>rba|id|name|set_name|stoichiometry')
        raise UserWarning('Invalid input file')
