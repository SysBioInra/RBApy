"""
Module defining UnknownProteins class.
"""

# python 2/3 compatibility
from __future__ import division, print_function

# global imports
import os.path
import pandas

# local imports
from rba.prerba.curation_data import CurationData

class UnknownProteins(object):
    """
    Class used to parse unknown protein file.
    """

    _helper_file = 'unknown_proteins.tsv'

    def __init__(self, input_dir='.'):
        """
        Constructor.

        Attributes:
            input_dir: path to folder containing helper files.
        """
        print('\nParsing unknown proteins'
              '\n------------------------')

        self._curated_data = CurationData(['SBML ID', 'UNIPROT GENE'])
        self._data_file = os.path.join(input_dir, self._helper_file)
        try:
            self._curated_data.read(self._data_file)
            print('Helper file found.')
        except IOError:
            print('Helper file not found.')
        self.data = {}
        for sbml, uniprot in self._curated_data.data.values:
            if pandas.isnull(uniprot):
                uniprot = None
            self.data[sbml] = uniprot

    def _print_missing_warning(self):
        print('WARNING: Several enzymatic proteins defined in your SBML '
              'could not be retrieved. '
              'Please read file {}, check data and specify all missing '
              'information. Execution will continue with default '
              'values (unknown proteins replaced by average enzymatic '
              'protein in the cytosol).'.format(self._data_file))

    def add(self, sbml_list):
        """
        Add identifiers that could not be retrieved to helper file.

        Args:
            sbml_list: list of sbml proteins that could not be retrieved in
                uniprot.
        """
        invalid_data = []
        data_to_cure = []
        for sbml in sbml_list:
            try:
                uniprot = self.data[sbml]
                if (uniprot != '' and not pandas.isnull(uniprot)
                    and not uniprot.startswith('average_protein_')):
                    invalid_data.append([sbml, self.data[sbml]])
            except KeyError:
                self.data[sbml] = None
                data_to_cure.append((sbml, None))

        if data_to_cure:
            self._curated_data.add(data_to_cure)
            self._curated_data.write(self._data_file)

        if invalid_data:
            print('WARNING: in {}, following lines do not map '
                  'to a valid uniprot gene, average or null protein:\n'
                  .format(self._data_file)
                  + '\n'.join(['\t'.join(d) for d in invalid_data]))

        if self._curated_data.has_missing_information():
            self._print_missing_warning()
