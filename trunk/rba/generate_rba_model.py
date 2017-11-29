"""Build RBA XML files from SBML."""

import sys

from rba import build_model


if __name__ == '__main__':
    # get parameter file
    if len(sys.argv) < 2:
        print('Please provide path to parameter file as script parameter.')
    else:
        try:
            model = build_model(sys.argv[1])
            model.write()
        except UserWarning as error:
            print(error.message)
            print('Interrupting Execution')
