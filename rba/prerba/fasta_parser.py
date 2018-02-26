"""Module defining FastaEntry and FastaParser classes."""

# python 2/3 compatibility
from __future__ import division, print_function, absolute_import

# global imports
from collections import namedtuple
from Bio import SeqIO

# class holding fasta entries
FastaEntry = namedtuple('FastaEntry',
                        'id name set_name stoichiometry sequence')


def parse_rba_fasta(input_file):
    """
    Parse rba-formatted fasta file.

    Parameters
    ----------
    input_file : str
        Path to fasta file.

    Returns
    -------
    list of FastaEntry
        Entries read

    """
    try:
        return [parse_entry(r) for r in SeqIO.parse(input_file, 'fasta')]

    except IOError:
        raise UserWarning('Please provide file ' + input_file)


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
