"""Build RBA XML files from SBML."""

import sys

from rba import ModelBuilder


if __name__ == '__main__':
    # get parameter file
    if len(sys.argv) < 2:
        print('Please provide path to parameter file as script parameter.')
    else:
        try:
            builder = ModelBuilder(sys.argv[1])
            builder.export_proteins('protein_summary.tsv')
            print('Building model...')
            model = builder.build_model()
            model.write()
            print('Done')
        except UserWarning as error:
            print(error)
            print('Interrupting Execution')
