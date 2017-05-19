
import os.path

class Macrocomponents(object):
    """
    Class used to parse/store target fluxes for macrocomponents.
    """
    def __init__(self, known_species, input_dir = '.'):
        """
        Constructor.
        :param input_dir: path to folder containing concentration file.
        :type input_dir: string
        """
        # constant values
        self._data_file = os.path.join(input_dir, 'macrocomponents.tsv')
        
        # read data
        self.target_flux = {}
        self._read_curated_data(known_species)
        
    def _read_curated_data(self, known_species):
        """
        Read file containing hand-curated information.
        """
        try:
            with open(self._data_file, 'rU') as input_stream:
                print('Found file with macrocomponent data. This file will be '
                      'used to setup macrocomponent production fluxes...')
                # skip header
                try:
                    next(input_stream)
                except StopIteration:
                    pass
                # read lines
                for line in input_stream:
                    line = line.rstrip('\n')
                    [id_, flux] = line.split('\t')
                    if id_ in known_species:
                        self.target_flux[id_] = float(flux)
                    else:
                        print('\nWARNING: in file ' + self._data_file + ': '
                              + id_ + ' does not correspond to a valid SBML '
                              'species. Line will be ignored.\n')
        except IOError:
            print('\nWARNING: Could not find file with macrocomponents '
                  'data... Run is going to continue without macrocomponent '
                  'synthesis. '
                  'Update helper file ' + self._data_file + ' with your own '
                  'values in future runs.\n')
            with open(self._data_file, 'w') as output_stream:
                # write header
                output_stream.write('\t'.join(['TARGET_METABOLITE',
                                               'TARGET_FLUX']) + '\n')
            
