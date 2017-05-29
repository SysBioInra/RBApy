
if __name__ == '__main__':
    # read data
    with open('unknown_proteins.tsv', 'r') as f:
        lines = f.readlines()

    # transform data
    new_lines = [lines[0]]
    for l in lines[1:]:
        sbml, uniprot = l.rstrip('\n').split('\t')
        # detect spontaneous reactions, remove RNAs
        # d -> diffusion
        # s -> spontaneous
        # RS0 -> RNA
        if any((sbml.startswith(prefix) for prefix in ['d', 's', 'RS0'])):
            uniprot = ''
        # replace unknown proteins with average proteins
        # e -> unknown excretion proteins
        elif sbml.startswith('e'):
            uniprot = 'average_protein_Cell_outer_membrane'
        # NoAssignment -> unknown proteins ???
        elif sbml == 'NoAssignment':
            uniprot = 'average_protein_Cytoplasm'
        new_lines.append('\t'.join([sbml, uniprot]) + '\n')

    # write transformed date to file
    with open('unknown_proteins.tsv', 'w') as f:
        f.write(''.join(new_lines))
