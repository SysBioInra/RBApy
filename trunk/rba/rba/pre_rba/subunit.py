
import os.path
import re
import pandas

class Subunit:
    """
    Class used to store curated subunit information.
    """
    def __init__(self, stoichiometry, gene_names, name, uniprot_note):
        """
        Constructor.

        :param stoichiometry: stoichiometry of protein within its complex.
        :param gene_names: Names of genes producing protein.
        :param name: Name of protein.
        :param uniprot_note: Uniprot note used to extract subunit information.
        :type stoichiometry: float
        :type gene_names: string
        :type name: string
        :type uniprot_note: string
        """
        self.stoichiometry = stoichiometry
        self.gene_names = gene_names
        self.name = name
        self.uniprot_note = uniprot_note

class SubunitData:
    """
    Class used to parse uniprot subunit field.
    """
    def __init__(self, uniprot_data, input_dir = '.'):
        """
        Constructor.

        :param uniprot_data: uniprot data as obtained by pandas.read_csv.
        :param input_dir: path to folder containing uniprot files.
        :type uniprot_data: pandas data
        :type input_dir: string
        """
        # constant values
        self._missing_string = ''
        self._missing_num = -1
        self._missing_tag = '[MISSING]'
        self._default_stoichiometry = 1
        self._subunit_file = os.path.join(input_dir, 'subunits.tsv')

        # load curated data
        self._missing_information = False
        self._curated_subunits = {}
        self._read_curated_data()

        # read uniprot data
        self.subunits = {}
        self._extract_uniprot_subunits(uniprot_data)

        # print warning if information was missing
        if self._missing_information:
            print('\nWARNING: Several uniprot notes were ambiguous. '
                  'Please read file '
                  + self._subunit_file
                  + ', check data and specify all information tagged as '
                  + self._missing_tag + '. Execution will continue with '
                  'default values (missing stoichiometry treated as '
                  + str(self._default_stoichiometry) + '). Make sure to update '
                  'helper file for next execution.\n')
            
    def _read_curated_data(self):
        """
        Read file containing hand-curated information.
        """
        try:
            with open(self._subunit_file, 'rU') as input_stream:
                print('Found file with subunit data. This file will be used '
                      'to solve ambiguous uniprot annotation...')
                # skip header
                next(input_stream)
                # read lines
                for line in input_stream:
                    [entry, stoichiometry,
                     gene_names, name, uniprot_note] = line.split('\t')
                    if stoichiometry == self._missing_tag:
                        self._missing_information = True
                        stoichiometry = self._missing_num
                    else:
                        stoichiometry = float(stoichiometry)
                    self._curated_subunits[entry] \
                        = Subunit(stoichiometry, gene_names, name, uniprot_note)
        except IOError:
            print('Could not find file with subunit data...')

    def _write_curated_data(self):
        """
        Write file containing hand-curated information.
        """
        with open(self._subunit_file, 'w') as output_stream:
            cols = ['ENTRY','STOICHIOMETRY', 'GENE NAMES',
                    'PROTEIN NAME', 'UNIPROT NOTE']
            output_stream.write('\t'.join(cols) + '\n')
            for entry in self._curated_subunits:
                subunit = self._curated_subunits[entry]
                if subunit.stoichiometry == self._missing_num:
                    sto = self._missing_tag
                else:
                    sto = str(subunit.stoichiometry)
                data = [entry, sto, subunit.gene_names,
                        subunit.name, subunit.uniprot_note]
                output_stream.write('\t'.join(data) + '\n')
            
    def _extract_uniprot_subunits(self, uniprot_data):
        """
        Extract subunit information from uniprot file, using curated data
        if available or writing missing data to file for future curation.
        """
        entry_data = uniprot_data['Entry']
        gene_data = uniprot_data['Gene names']
        name_data = uniprot_data['Protein names']
        subunit_data = uniprot_data['Subunit structure [CC]']

        subunit_reader = re.compile(r'([a-zA-Z]+)mer[^a-z]')

        curated_data_added = False
        for entry, gene_names, name, subunit_field \
            in zip(entry_data, gene_data, name_data, subunit_data):
            # if entry is in curated file, simply add it
            if self._curated_subunits.has_key(entry):
                stoichiometry = self._curated_subunits[entry].stoichiometry
                if stoichiometry == self._missing_num:
                    self.subunits[entry] = self._default_stoichiometry
                else:
                    self.subunits[entry] = stoichiometry
                continue
            
            # if field is empty, stochiometry is one
            if pandas.isnull(subunit_field):
                self.subunits[entry] = 1
                continue
            # else try to deduce stoichiometry from subunit field
            prefixes = subunit_reader.findall(subunit_field)
            # if there is only one word of the form [prefix]mer,
            if len(prefixes) == 1:
                prefix = prefixes[0].upper()
                # we automatically accept the following:
                stoichiometry = 0
                if prefix == 'MONO':
                    stoichiometry = 1
                elif prefix == 'HETERODI':
                    stoichiometry = 1
                elif prefix == 'HOMODI':
                    stoichiometry = 2
                elif prefix == 'HOMOTRI':
                    stoichiometry = 3
                elif prefix == 'HOMOTETRA':
                    stoichiometry = 4
                elif prefix == 'HOMOPENTA':
                    stoichiometry = 5
                elif prefix == 'HOMOHEXA':
                    stoichiometry = 6
                elif prefix == 'HEPTA':
                    stoichiometry = 7
                elif prefix == 'HOMOOCTA':
                    stoichiometry = 8
                elif prefix == 'HOMODECA':
                    stoichiometry = 10
                elif prefix == 'HOMODODECA':
                    stoichiometry = 12
                if stoichiometry > 0:
                    self.subunits[entry] = stoichiometry
                    continue

            # if we get here, field was ambiguous:
            # tag all missing info and write to curated file
            curated_data_added = True
            self._missing_information = True
            self._curated_subunits[entry] \
                = Subunit(self._missing_num, gene_names, name, subunit_field)
            # use default stoichiometry for the time being
            self.subunits[entry] = self._default_stoichiometry

        if curated_data_added:
            self._write_curated_data()
