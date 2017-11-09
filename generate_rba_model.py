"""Build RBA XML files from SBML."""

import sys

from rba import PreRba


if __name__ == '__main__':
    # get parameter file
    if len(sys.argv) < 2:
        print('Please provide path to parameter file as script parameter.')
    else:
        try:
            pre_rba = PreRba(sys.argv[1])
            pre_rba.model.write_files()
        except UserWarning as error:
            print(error.message)
            print('Interrupting Execution')
