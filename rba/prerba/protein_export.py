"""Module exporting protein to a flat file."""

# python 2/3 compatibility
from __future__ import division, print_function, absolute_import


def export_proteins(filename, protein_data):
    with open(filename, 'w') as output_file:
        output_file.write(header() + '\n')
        output_file.write('\n'.join(format_protein(p) for p in protein_data)
                          + '\n')


def header():
    return '\t'.join(['IDENTIFIER', 'SEQUENCE', 'COFACTORS', 'LOCATION',
                      'STOICHIOMETRY'])


def format_protein(protein):
    cofactors = '; '.join(format_cofactor(c) for c in protein.cofactors)
    return '\t'.join([protein.id, protein.sequence, cofactors,
                      protein.location, str(protein.stoichiometry)])


def format_cofactor(cofactor):
    return '{} {} ({})'.format(cofactor.stoichiometry, cofactor.chebi,
                               cofactor.name)
