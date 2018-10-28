"""Build RBA XML files from SBML."""

import sys

from rba import RbaModel


if __name__ == '__main__':
    # get parameter file
    if len(sys.argv) < 2:
        print('Please provide path to parameter file as script parameter.')
    else:
        try:
            model = RbaModel.from_data(sys.argv[1])
            model.write()
            print('Done')
        except UserWarning as error:
            print(error)
            print('Interrupting Execution')
