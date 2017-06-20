
"""
Module defining Cofactor, CofactorParser and CofactorData classes.
"""

# python 2/3 compatibility
from __future__ import absolute_import, division, print_function

# global imports
import os.path
import re
from collections import namedtuple
import pandas

# Class storing cofactor information.
Cofactor = namedtuple('Cofactor', 'chebi name stoichiometry uniprot_note')

class CofactorParser(object):
    """
    Class parsing Cofactor uniprot field.
    """

    def __init__(self):
        """
        Constructor.
        """
        self._name_reader = re.compile(r'Name=([^;]+); Xref=ChEBI:([^;]+);')
        self._note_reader = re.compile(r'Note=(.*)')
        self._stoichiometry_reader = re.compile(r'Binds ([\w]+)')

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

class CofactorData(object):
    """
    Class parsing uniprot cofactor field.

    Attributes:
        cofactors: dict mapping uniprot entry with list of Cofactor objects.
    """

    def __init__(self, uniprot_data, input_dir='.'):
        """
        Constructor from pandas dataframe.

        Args:
            uniprot_data: pandas dataframe containing uniprot data.
            input_dir: path to folder containing uniprot file.
        """
        # constant values
        self._missing_tag = '[MISSING]'
        self._default_stoichiometry = 1

        # read curated cofactors
        curation_file = os.path.join(input_dir, 'cofactors.tsv')
        try:
            with open(curation_file, 'rU') as input_stream:
                print('Found file with cofactor data. This file will be used'
                      ' to solve ambiguous uniprot annotation...')
                # skip header
                next(input_stream)
                curated_data, missing_information \
                    = self._read_curated_data(input_stream)
        except IOError:
            print('Could not find file with cofactor data...')
            curated_data = {}
            missing_information = False

        # read uniprot data
        self.cofactors, cofactors_to_cure \
            = self._extract_data(uniprot_data, curated_data)
        self.cofactors.update(cofactors_to_cure)
        self.cofactors = self._fill_missing_information(self.cofactors)

        # write curated data to file (if necessary)
        if cofactors_to_cure:
            curated_data.update(cofactors_to_cure)
            with open(curation_file, 'w') as output_stream:
                self._write_curated_data(output_stream, curated_data)

        # print warning if data was missing
        if missing_information or cofactors_to_cure:
            print('\nWARNING: Several uniprot notes were ambiguous. '
                  'Please read file {}, check data and specify all '
                  'information tagged as {}. Add lines if necessary but avoid '
                  'removing lines. Execution will continue with default '
                  'values (missing stoichiometry treated as {}). '
                  'Make sure to update helper file for next execution.\n'
                  .format(curation_file, self._missing_tag,
                          self._default_stoichiometry))

    def _fill_missing_information(self, cofactors):
        """
        Assign default values to missing fields.
        """
        new_cofactors = {}
        for entry, cofactor_list in cofactors.items():
            for cofactor in cofactor_list:
                # ignore cofactors with unknown chebi or zero stoichiometry
                if cofactor.chebi is not None and cofactor.stoichiometry != 0:
                    # replace missing stoichiometry with default stoichiometry
                    if cofactor.stoichiometry is None:
                        new_cofactor = cofactor._replace(stoichiometry=self._default_stoichiometry)
                    else:
                        new_cofactor = cofactor
                    new_cofactors.setdefault(entry, []).append(new_cofactor)
        return new_cofactors

    def _read_curated_data(self, input_stream):
        """
        Read file containing hand-curated cofactor information.
        """
        missing_information = False
        curated_data = {}
        # read lines
        for line in input_stream:
            [entry, chebi, name,
             stoichiometry, uniprot_note] = line.split('\t')
            if chebi == self._missing_tag:
                missing_information = True
                chebi = None
            if stoichiometry == self._missing_tag:
                missing_information = True
                stoichiometry = None
            else:
                stoichiometry = float(stoichiometry)
            new_cofactor = Cofactor(chebi, name, stoichiometry, uniprot_note)
            curated_data.setdefault(entry, []).append(new_cofactor)
        return curated_data, missing_information

    def _write_curated_data(self, output_stream, curated_data):
        """
        Write file containing hand-curated cofactor information.
        """
        cols = ['ENTRY', 'CHEBI', 'NAME',
                'STOICHIOMETRY', 'UNIPROT ANNOTATION']
        output_stream.write('\t'.join(cols) + '\n')
        for entry, cofactors in curated_data.items():
            for cofactor in cofactors:
                if cofactor.chebi is None:
                    name = chebi = self._missing_tag
                else:
                    name = cofactor.name
                    chebi = cofactor.chebi
                if cofactor.stoichiometry is None:
                    sto = self._missing_tag
                else:
                    sto = str(cofactor.stoichiometry)
                data = [entry, chebi, name, sto, cofactor.uniprot_note]
                output_stream.write('\t'.join(data) + '\n')

    @staticmethod
    def _extract_data(uniprot_data, curated_data):
        """
        Extract cofactor information from uniprot file.
        """
        parser = CofactorParser()
        entry_data = uniprot_data['Entry']
        cofactor_data = uniprot_data['Cofactor']
        cofactors = {}
        cofactors_to_cure = {}
        for entry, cofactor_field in zip(entry_data, cofactor_data):
            # if entry belongs to curated set, simply add curated data
            curated_cofactors = curated_data.get(entry, None)
            if curated_cofactors is not None:
                cofactors[entry] = curated_cofactors
            else:
                if pandas.isnull(cofactor_field):
                    continue
                cofactors[entry], to_cure = parser.parse(cofactor_field)
                if to_cure:
                    cofactors_to_cure[entry] = to_cure
        return cofactors, cofactors_to_cure
