
import os.path
import re
import pandas

class Cofactor:
    """
    Class used to store cofactor information.
    """
    def __init__(self, chebi, name, stoichiometry, uniprot_note = ''):
        """
        Constructor.

        :param chebi: ChEBI identifier of cofactor.
        :param name: Name of cofactor.
        :param stoichiometry: Stoichiometry of cofactor.
        :param uniprot_note: Uniprot note used to extract cofator information.
        :type chebi: string
        :type name: string
        :type stoichiometry: float
        :type uniprot_note: string
        """
        self.chebi = chebi
        self.name = name
        self.stoichiometry = stoichiometry
        self.uniprot_note = uniprot_note

class CofactorData:
    """
    Class used to parse uniprot cofactor field.
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
        self._cofactor_file = os.path.join(input_dir, 'cofactors.tsv')
        
        # read curated cofactors
        self._missing_information = False
        self._curated_cofactors = {}
        self._read_curated_data()

        # read uniprot data
        self.cofactors = {}
        self._extract_uniprot_cofactors(uniprot_data) \
        
        # raise warning if data was missing
        if self._missing_information:
            print('\nWARNING: Several uniprot notes were ambiguous. '
                  'Please read file ' + self._cofactor_file
                  + ', check data and specify all information tagged as '
                  + self._missing_tag + '. Add lines if necessary but avoid '
                  'removing lines. Execution will continue with default '
                  'values (missing stoichiometry treated as '
                  + str(self._default_stoichiometry) + '). Make sure to update '
                  'helper file for next execution.\n')

    def _add_cofactor(self, entry, chebi, name, stoichiometry):
        """
        Add cofactor to dictionary of read cofactors.
        """
        # ignore cofactors with unknown names
        if (chebi == self._missing_string): return
        # replace missing stoichiometry with default stoichiometry
        if (stoichiometry == self._missing_num):
            stoichiometry = self._default_stoichiometry
        # ignore cofactors with zero stoichiometry
        if (stoichiometry <= 0): return
            
        if entry not in self.cofactors:
            self.cofactors[entry] = []
        new_cofactor = Cofactor(chebi, name, stoichiometry);
        self.cofactors[entry].append(new_cofactor)

    def _add_curated_cofactor(self, entry, chebi, name, stoichiometry, uniprot_note):
        """
        Add curated cofactor to dictionary of curated cofactors.
        """
        if entry not in self._curated_cofactors:
            self._curated_cofactors[entry] = []
        new_cofactor = Cofactor(chebi, name, stoichiometry, uniprot_note);
        self._curated_cofactors[entry].append(new_cofactor)

    def _read_curated_data(self):
        """
        Read file containing hand-curated cofactor information.
        """
        try:
            with open(self._cofactor_file, 'rU') as input_stream:
                print('Found file with cofactor data. This file will be used'
                      ' to solve ambiguous uniprot annotation...')
                # skip header
                next(input_stream)
                # read lines
                for line in input_stream:
                    [entry, chebi, name, stoichiometry, uniprot_note] \
                        = line.split('\t')
                    if chebi == self._missing_tag:
                        self._missing_information = True
                        chebi = self._missing_string
                    if stoichiometry == self._missing_tag:
                        self._missing_information = True
                        stoichiometry = self._missing_num
                    else:
                        stoichiometry = float(stoichiometry)
                    self._add_curated_cofactor(entry, chebi, name,
                                               stoichiometry, uniprot_note)
        except IOError:
            print('Could not find file with cofactor data...')

    def _write_curated_data(self):
        """
        Write file containing hand-curated cofactor information.
        """
        with open(self._cofactor_file, 'w') as output_stream:
            cols = ['ENTRY', 'CHEBI', 'NAME',
                    'STOICHIOMETRY', 'UNIPROT ANNOTATION']
            output_stream.write('\t'.join(cols) + '\n')
            for entry in self._curated_cofactors:
                for cofactor in self._curated_cofactors[entry]:
                    if cofactor.name == self._missing_string:
                        name = chebi = self._missing_tag
                    else:
                        name = cofactor.name
                        chebi = cofactor.chebi
                    if cofactor.stoichiometry == self._missing_num:
                        sto = self._missing_tag
                    else:
                        sto = str(cofactor.stoichiometry)
                    data = [entry, chebi, name, sto, cofactor.uniprot_note]
                    output_stream.write('\t'.join(data) + '\n')
        
    def _extract_uniprot_cofactors(self, uniprot_data):
        """
        Extract cofactor information from uniprot file, using curated data
        if available or writing missing data to file for future curation.
        """
        entry_data = uniprot_data['Entry']
        cofactor_data = uniprot_data['Cofactor']

        cofactor_reader = re.compile(r'Name=([^;]+); Xref=ChEBI:([^;]+);')
        note_reader = re.compile(r'Note=(.*)')
        stoichiometry_reader = re.compile(r'Binds ([\w]+)')

        curated_data_added = False
        for entry, cofactor_field in zip(entry_data, cofactor_data):
            # check whether entry belongs to curated set
            if entry in self._curated_cofactors:
                self.cofactors[entry] = []
                for c in self._curated_cofactors[entry]:
                    self._add_cofactor(entry, c.chebi, c.name, c.stoichiometry)
                continue
            
            # split different cofactors
            if pandas.isnull(cofactor_field): continue
            cofactor_notes = cofactor_field.split('COFACTOR:')[1:]
            if len(cofactor_notes) == 0: continue

            # read each cofactor note
            for note in cofactor_notes:
                full_name = cofactor_reader.findall(note)
                stoichiometry = stoichiometry_reader.findall(note)
                # infer stoichiometry
                if len(stoichiometry) == 0:
                    stoichiometry = 1
                elif len(stoichiometry) == 1:
                    try:
                        stoichiometry = int(stoichiometry[0])
                    except ValueError:
                        # number is probably given as a string, but we just
                        # tag it as missing
                        stoichiometry = self._missing_num
                else: # ambiguous case
                    stoichiometry = self._missing_num
                    
                # if there is one cofactor and stoichiometry isn't missing
                # there is no ambiguity
                if (len(full_name) == 1) and (stoichiometry != self._missing_num):
                    [name, chebi] = full_name[0]
                    self._add_cofactor(entry, chebi, name, stoichiometry)
                # else result is ambiguous: we need help from the user:
                # tag all missing info and write to curated file
                else:
                    self._missing_information = True
                    curated_data_added = True
                    # if we could not find a name, indicate chebi and name
                    # as missing
                    if len(full_name) == 0:
                        full_name.append([self._missing_string, \
                                          self._missing_string])
                    # extract subnote if possible and display instead
                    # of whole note (if extraction of subnote is ambiguous,
                    # we display the whole note).
                    subnotes = note_reader.findall(note)
                    if len(subnotes) == 1:
                        note = subnotes[0]
                    # if there are several names, assume stoichiometry
                    # is number found earlier for first element of the list
                    # and 0 for the rest
                    for name, chebi in full_name:
                        # add cofactor (missing values will be replaced
                        # by default values)
                        self._add_cofactor(entry, chebi, name, stoichiometry)
                        # add curation note
                        self._add_curated_cofactor(entry, chebi, name, \
                                                   stoichiometry, note)
                        stoichiometry = 0

        # write curated data to file (if necessary)
        if curated_data_added: self._write_curated_data()
