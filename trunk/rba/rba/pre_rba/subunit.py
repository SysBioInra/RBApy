"""
Module defining Subunit, SubunitParser and SubunitData classes.
"""

# python 2/3 compatibility
from __future__ import division, print_function

# global imports
import os.path
import re
from collections import namedtuple
import pandas

# local imports
from rba.pre_rba.curation_data import CurationData

# Class used to store curated subunit information.
Subunit = namedtuple('Subunit', 'stoichiometry gene_names name uniprot_note')

class SubunitData(object):
    """
    Class used to parse uniprot subunit field.
    """
    default_stoichiometry = 1
    _helper_file = 'subunits.tsv'

    def __init__(self, uniprot_data, input_dir='.'):
        """
        Constructor from pandas.DataFrame.

        Args:
             uniprot_data: uniprot data as obtained by pandas.read_csv.
             input_dir: path to folder containing uniprot files.
        """
        print('\nParsing subunits'
              '\n----------------')
        # load curated data
        curation_data = CurationData(['ENTRY', 'STOICHIOMETRY', 'GENE NAMES',
                                      'PROTEIN NAME', 'UNIPROT NOTE'])
        curation_file = os.path.join(input_dir, self._helper_file)
        try:
            curation_data.read(curation_file)
            print('Helper file found.')
        except IOError:
            print('Helper file not found.')
        curated_subunits = {row[0]: Subunit(*row[1:])
                            for row in curation_data.data.values}

        # read uniprot data
        # store new data that needs curation (if applicable)
        self.subunits, data_to_cure \
            = self._extract_data(uniprot_data, curated_subunits)
        if data_to_cure:
            curation_data.add(data_to_cure)
            curation_data.write(curation_file)

        # print warning if information was missing
        if curation_data.has_missing_information():
            print('WARNING: Several uniprot notes were ambiguous. '
                  'Please read file {}, check data and specify all '
                  'information tagged as {}. Execution will continue with '
                  'default values (missing stoichiometry treated as {}).'
                  .format(curation_file, CurationData.missing_tag,
                          self.default_stoichiometry))

    def _extract_data(self, uniprot_data, curated_subunits):
        """
        Extract subunit information from uniprot file.
        """
        entry_data = uniprot_data['Entry']
        gene_data = uniprot_data['Gene names']
        name_data = uniprot_data['Protein names']
        subunit_data = uniprot_data['Subunit structure [CC]']
        parser = SubunitParser()

        subunits = {}
        data_to_cure = []
        for entry, gene_names, name, subunit_field \
            in zip(entry_data, gene_data, name_data, subunit_data):
            # if entry is in curated file, add it, otherwise parse field
            if entry in curated_subunits:
                stoichiometry = curated_subunits[entry].stoichiometry
                if pandas.isnull(stoichiometry):
                    subunits[entry] = self.default_stoichiometry
                else:
                    subunits[entry] = stoichiometry
            else:
                stoichiometry = parser.parse(subunit_field)
                if stoichiometry is not None:
                    subunits[entry] = stoichiometry
                else:
                    # field was ambiguous:
                    # use default stoichiometry, add to data to cure
                    data_to_cure.append((entry, None, gene_names,
                                         name, subunit_field))
                    subunits[entry] = self.default_stoichiometry
        return subunits, data_to_cure

class SubunitParser(object):
    """
    Class parsing 'Subunit' uniprot field.

    Attributes:
        prefix_rule: dictionary determining rule used to infer stoichiometry.
            Keys are all caps prefixes preceding 'mer' in words found
            in uniprot field,
            values are stoichiometries associated with them. For example,
            prefix_rule[MONO] = 1.
    """

    prefix_rule = {'MONO': 1, 'HETERODI': 1, 'HOMODI': 2, 'HOMOTRI': 3,
                   'HOMOTETRA': 4, 'HOMOPENTA': 5, 'HOMOHEXA': 6,
                   'HEPTA': 7, 'HOMOOCTA': 8, 'HOMODECA': 10, 'HOMODODECA': 12}
    _subunit_reader = re.compile(r'([a-zA-Z]+)mer[^a-z]')

    def parse(self, field):
        """
        Parse uniprot field.

        Args:
            field: field to parse.

        Returns:
            Stoichiometry parsed (None if field was ambiguous).
        """
        # if field is empty, stoichiometry is one
        if pandas.isnull(field): return 1

        prefixes = self._subunit_reader.findall(field)
        # if there is only one word of the form [prefix]mer,
        if len(prefixes) == 1:
            prefix = prefixes[0].upper()
            # check for prefix in prefix rule
            return self.prefix_rule.get(prefix, None)
        else:
            return None
