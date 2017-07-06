
"""
Module defining Cofactor, CofactorParser and CofactorData classes.
"""

# python 2/3 compatibility
from __future__ import division, print_function

# global imports
import os.path
import re
from collections import namedtuple
import pandas

# local imports
from rba.prerba.curation_data import CurationData

# Class storing cofactor information.
Cofactor = namedtuple('Cofactor', 'chebi name stoichiometry uniprot_note')

class CofactorData(object):
    """
    Class parsing uniprot cofactor field.

    Attributes:
        default_stoichiometry: stoichiometry assigned by default.
        cofactors: dict mapping uniprot entry with list of Cofactor objects.
    """

    default_stoichiometry = 1
    _helper_file = 'cofactors.tsv'

    def __init__(self, uniprot_data, input_dir='.'):
        """
        Constructor from pandas dataframe.

        Args:
            uniprot_data: pandas dataframe containing uniprot data.
            input_dir: path to folder containing uniprot file.
        """
        print('\nParsing cofactors'
              '\n-----------------')
        # import curation data
        curation_data = CurationData(['ENTRY', 'CHEBI', 'NAME',
                                      'STOICHIOMETRY', 'UNIPROT ANNOTATION'])
        curation_file = os.path.join(input_dir, self._helper_file)
        try:
            curation_data.read(curation_file)
            print('Helper file found.')
        except IOError:
            print('Helper file not found.')
        curated_cofactors = {}
        for row in curation_data.data.values:
            curated_cofactors.setdefault(row[0], []).append(Cofactor(*row[1:]))

        # read uniprot data
        # if new data to cure is found, add it to curation file
        self.cofactors, data_to_cure \
            = self._extract_data(uniprot_data, curated_cofactors)
        if data_to_cure:
            curation_data.add(data_to_cure)
            curation_data.write(curation_file)

        self.cofactors = self._fill_missing_information(self.cofactors)

        # print warning if data needs curation
        if curation_data.has_missing_information():
            print('WARNING: Several uniprot notes were ambiguous. '
                  'Please read file {}, check data and specify all '
                  'information tagged as {}. Add lines if necessary but avoid '
                  'removing lines. Execution will continue with default '
                  'values (missing stoichiometry treated as {}).'
                  .format(curation_file, CurationData.missing_tag,
                          self.default_stoichiometry))

    def _fill_missing_information(self, cofactors):
        """
        Assign default values to missing fields.
        """
        new_cofactors = {}
        for entry, cofactor_list in cofactors.items():
            for cofactor in cofactor_list:
                # ignore cofactors with unknown chebi or zero stoichiometry
                if (not pandas.isnull(cofactor.chebi)
                    and cofactor.stoichiometry != 0):
                    # replace missing stoichiometry with default stoichiometry
                    if pandas.isnull(cofactor.stoichiometry):
                        new_cofactor = cofactor._replace(stoichiometry=self.default_stoichiometry)
                    else:
                        new_cofactor = cofactor
                    new_cofactors.setdefault(entry, []).append(new_cofactor)
        return new_cofactors

    @staticmethod
    def _extract_data(uniprot_data, curated_cofactors):
        """
        Extract cofactor information from uniprot file.
        """
        parser = CofactorParser()
        entry_data = uniprot_data['Entry']
        cofactor_data = uniprot_data['Cofactor']
        cofactors = {}
        data_to_cure = []
        for entry, cofactor_field in zip(entry_data, cofactor_data):
            # if entry belongs to curated set, simply add curated data
            curated = curated_cofactors.get(entry, None)
            if curated is not None:
                cofactors[entry] = curated
            else:
                if pandas.isnull(cofactor_field):
                    continue
                new_cofactors, to_cure = parser.parse(cofactor_field)
                cofactors[entry] = new_cofactors + to_cure
                data_to_cure += [(entry,) + c for c in to_cure]
        return cofactors, data_to_cure

class CofactorParser(object):
    """
    Class parsing Cofactor uniprot field.
    """

    _name_reader = re.compile(r'Name=([^;]+); Xref=ChEBI:([^;]+);')
    _note_reader = re.compile(r'Note=(.*)')
    _stoichiometry_reader = re.compile(r'Binds ([\w]+)')

    def parse(self, field):
        """
        Parse uniprot field.

        Args:
            field: string containing uniprot cofactor information.

        Returns:
            Two elements. First element is a list of Cofactor objects
            containing information retrieved, where info was unambiguous.
            Second element is a list of Cofactor objects where some info was
            ambiguous. If some information could not be retrieved,
            its field is set to None
        """
        cofactor_notes = field.split('COFACTOR:')[1:]
        cofactors = []
        cofactors_to_cure = []
        for note in cofactor_notes:
            # read name(s) and chebi identifier(s) of cofactor(s)
            # if no name was found, indicate chebi and name as missing
            full_name = self._name_reader.findall(note)
            if not full_name:
                full_name.append([None, None])
            # extract subnote if possible
            subnotes = self._note_reader.findall(note)
            if len(subnotes) == 1:
                subnote = subnotes[0]
            else:
                subnote = note
            # infer stoichiometry:
            #  - nothing read: stoichiometry is implicitly 1
            #  - one value read: use value if can be cast to integer, else
            #    tag as missing information.
            #  - multiple values read: tag as missing information.
            stoichiometry = self._stoichiometry_reader.findall(note)
            if not stoichiometry:
                stoichiometry = 1
            elif len(stoichiometry) == 1:
                try:
                    stoichiometry = int(stoichiometry[0])
                except ValueError:
                    stoichiometry = None
            else:
                stoichiometry = None
            # if there are several names, assume stoichiometry
            # is number found earlier for first element of the list
            # and 0 for the rest
            is_ambiguous = (stoichiometry is None
                            or len(full_name) > 1
                            or full_name[0][0] is None)
            for name, chebi in full_name:
                new_cofactor = Cofactor(chebi, name, stoichiometry, subnote)
                if is_ambiguous:
                    cofactors_to_cure.append(new_cofactor)
                else:
                    cofactors.append(new_cofactor)
                stoichiometry = 0
        return cofactors, cofactors_to_cure
