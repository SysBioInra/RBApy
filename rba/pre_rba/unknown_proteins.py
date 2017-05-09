
import os.path

class UnknownProteins:
    """
    Class used to parse uniprot cofactor field.
    """
    def __init__(self, input_dir = '.'):
        """
        Constructor.

        :param input_dir: path to folder containing helper files.
        """
        # constant values
        self._missing_tag = '[MISSING]'
        self._data_file = os.path.join(input_dir, 'unknown_proteins.tsv')
        
        # read curated cofactors
        self._missing_information = False
        self.data = {}
        self._read_curated_data()

    def _print_missing_warning(self):
        print('\nWARNING: Several enzymatic proteins defined in your SBML '
              'could not be retrieved. '
              'Please read file ' + self._data_file +
              ', check data and specify all information tagged as ' +
              self._missing_tag + '. Execution will continue with default '
              'values (unknown proteins replaced by average enzymatic '
              'protein in the cytosol). Make sure to update '
              'helper file for next execution.')

    def _read_curated_data(self):
        """
        Read file containing hand-curated cofactor information.
        """
        try:
            input_stream = open(self._data_file, 'r')
            print('Found file with unknown protein data. This file will be '
                  'used to match ambiguous enzymatic sbml annotations...')
            # skip header
            next(input_stream)
            # read lines
            for line in input_stream:
                [sbml, uniprot] = line.rstrip('\n').split('\t')
                if uniprot == self._missing_tag:
                    self._missing_information = True
                self.data[sbml] = uniprot
            input_stream.close()
        except IOError:
            print 'Could not find file with unknown protein data...'

    def _write_curated_data(self):
        """
        Write file containing hand-curated cofactor information.
        """
        output_stream = open(self._data_file, 'w')
        output_stream.write('\t'.join(['SBML ID', 'UNIPROT GENE']) + '\n')
        for sbml, uniprot in self.data.iteritems():
            output_stream.write('\t'.join([sbml, uniprot]) + '\n')
        output_stream.close()
        
    def add(self, sbml_list):
        """
        Add identifiers that could not be retrieved to hepler file.
        """
        missing_information = False
        invalid_data = []
        for sbml in sbml_list:
            try:
                uniprot = self.data[sbml]
                if uniprot != '' and uniprot != self._missing_tag \
                   and not(uniprot.startswith('average_protein_')):
                    invalid_data.append([sbml, self.data[sbml]])
            except KeyError:
                self.data[sbml] = self._missing_tag
                missing_information = True
        if len(invalid_data) > 0:
            print('\nWARNING: unknown_proteins.tsv, following lines do not map '
                  'to a valid uniprot gene, average or null protein:\n'
                  + '\n'.join(map('\t'.join, invalid_data)))
        if missing_information:
            self._write_curated_data()
            self._print_missing_warning()
        
