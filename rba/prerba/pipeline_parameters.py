"""Module defining PipelineParameters class."""

# python 2/3 compatibility
from __future__ import division, print_function, absolute_import


class PipelineParameters(object):
    """
    Class storing pipeline parameters.

    Attributes
    ----------
    obligatory_tags : list of str
        tags that a parameter file must contain.
    optional_tags : list of str
        tags that a parameter file may contain.
    parameters: dict
        Dictonary mapping parameter tags with their value. Optional tags
        may be omitted.

    """

    def __init__(self, parameter_file):
        """
        Build object from file path.

        Parameters
        ----------
        parameter_file: str
            path to parameter file.

        """
        self.obligatory_tags = ['INPUT_DIR', 'OUTPUT_DIR', 'SBML_FILE',
                                'ORGANISM_ID']
        self.optional_tags = ['EXTERNAL_COMPARTMENTS','INTERFACE_COMPARTMENTS']
        self.parameters = {}

        try:
            with open(parameter_file, 'rU') as input_stream:
                # parse file
                for line in input_stream:
                    line = line.strip()
                    if line == '' or line.startswith('#'):
                        continue
                    try:
                        tag, value = [word.strip() for word in line.split('=')]
                    except ValueError:
                        print ('Invalid parameter format:\n' + line)
                        raise UserWarning('Invalid parameter file.')
                    if (tag in self.obligatory_tags
                            or tag in self.optional_tags):
                        self.parameters[tag] = value
                    else:
                        print('Warning: ignoring unknown parameter '
                              + tag + '.')
        except IOError:
            print('Could not find parameter file ' + parameter_file + '.')
            raise UserWarning('Invalid parameter file.')

        # check that all obligatory tags have been found
        missing_parameters = [tag for tag in self.obligatory_tags
                              if tag not in self.parameters.keys()]
        if len(missing_parameters) > 0:
            print('Following parameters are missing: '
                  + ', '.join(missing_parameters))
            raise UserWarning('Invalid parameter file.')
