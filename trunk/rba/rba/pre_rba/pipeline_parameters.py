
class PipelineParameters:
    """
    Class storing pipeline parameters.
    """
    def __init__(self, parameter_file):
        """
        Constructor.

        :param parameter_file: path to parameter file.
        :type parameter_file: string
        """
        self.obligatory_tags = ['INPUT_DIR', 'OUTPUT_DIR', 'SBML_FILE',
                                'ORGANISM_ID']
        self.optional_tags = ['EXTERNAL_COMPARTMENTS']
        self.parameters = {}
        
        try:
            f = open(parameter_file)
        except IOError:
            print('Could not find parameter file ' + parameter_file + '.')
            raise UserWarning('Invalid parameter file.')

        # parse file
        for line in f:
            line = line.strip()
            if line == '' or line.startswith('#'): continue
            try:
                tag, value = map(str.strip, line.split('='))
            except ValueError:
                print ('Invalid parameter format:\n' + line)
                raise UserWarning('Invalid parameter file.')
            if (tag in self.obligatory_tags) \
               or (tag in self.optional_tags):
                self.parameters[tag] = value
            else:
                print('Warning: ignoring unknown parameter ' + tag + '.')
            
        # check that all obligatory tags have been found
        missing_parameters = [tag for tag in self.obligatory_tags
                              if tag not in self.parameters.keys()]
        if len(missing_parameters)>0:
            print('Following parameters are missing: '
                  + ', '.join(missing_parameters))
            raise UserWarning('Invalid parameter file.')
