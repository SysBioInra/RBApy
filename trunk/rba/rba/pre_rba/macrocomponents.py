"""
Module defining Macrocomponents class.
"""

# python 2/3 compatibility
from __future__ import division, print_function

# global imports
import os.path

# local imports
from rba.pre_rba.curation_data import CurationData

class Macrocomponents(object):
    """
    Class used to parse/store target fluxes for macrocomponents.

    Attributes:
        target_flux: dict mapping sbml species id to target production flux.
    """

    _helper_file = 'macrocomponents.tsv'

    def __init__(self, known_species, input_dir='.'):
        """
        Constructor.

        Args:
            known_species: list of sbml species ids.
            input_dir: path to folder containing concentration file.
        """
        print('\nParsing macrocomponent information'
              '\n----------------------------------')

        # read curation data
        curation_file = os.path.join(input_dir, self._helper_file)
        curation_data = CurationData(['TARGET_METABOLITE', 'TARGET_FLUX'])
        try:
            curation_data.read(curation_file)
            print('Helper file found.')
        except IOError:
            print('WARNING: Could not find file with macrocomponents '
                  'data... Run is going to continue without macrocomponent '
                  'synthesis. Use helper file {} in future runs.'
                  .format(curation_file))
            curation_data.write(curation_file)

        # read data
        self.target_flux = {}
        for metabolite, flux in curation_data.data.values:
            if metabolite in known_species:
                self.target_flux[metabolite] = float(flux)
            else:
                print('WARNING: in file {}: {} is not a valid SBML species. '
                      'Line will be ignored.'.format(curation_file, metabolite))
