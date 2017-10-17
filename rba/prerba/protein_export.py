"""Module exporting protein to a flat file."""


def export_proteins(filename, protein_data):
    with open(filename, 'w') as output_file:
        output_file.write(header() + '\n')
        output_file.write('\n'.join(format_protein(id_, p)
                                    for id_, p in protein_data.items())
                          + '\n')


def header():
    return '\t'.join(['IDENTIFIER', 'SEQUENCE', 'COFACTORS', 'LOCATION',
                      'STOICHIOMETRY'])


def format_protein(identifier, protein):
    cofactors = '; '.join(format_cofactor(c) for c in protein.cofactors)
    return '\t'.join([identifier, protein.sequence, cofactors,
                      protein.location, str(protein.stoichiometry)])


def format_cofactor(cofactor):
    return '{} {} ({})'.format(cofactor.stoichiometry, cofactor.chebi,
                               cofactor.name)
