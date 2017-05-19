
import os.path
import re
import pandas

class Location:
    """
    Class used to store curated location information.
    """
    def __init__(self, location, gene_names, name):
        """
        Constructor.

        :param location: usual location of protein.
        :param gene_names: Names of genes producing protein.
        :param name: Name of protein.
        :type location: string
        :type gene_names: string
        :type name: string
        """
        self.location = location
        self.gene_names = gene_names
        self.name = name

class LocationData:
    """
    Class used to parse uniprot location field.
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
        self.default_location = 'Cytoplasm'
        self._missing_string = ''
        self._missing_num = -1
        self._missing_tag = '[MISSING]'
        self._mapping_file = os.path.join(input_dir, 'location_map.tsv')
        self._curation_file = os.path.join(input_dir, 'locations.tsv')

        # load curated data
        self._missing_mapping = False
        self._location_map = {}
        self._read_mapping_data()
        self._missing_information = False
        self._curated_locations = {}
        self._read_curated_data()

        # read uniprot data
        self.locations = {}
        self._extract_uniprot_data(uniprot_data)

        # print warning if information is missing
        if self._missing_mapping:
            print('\nWARNING: Some uniprot locations had no user-defined '
                  'counterpart. '
                  + 'RBA compartments will reflect uniprot annotation.'
                  + 'If you wish to change the compartment names, please '
                  + 'read file' + self._mapping_file
                  + ' and specify all information tagged as '
                  + self._missing_tag + '.\n')

        # print warning if information is missing
        if self._missing_information:
            print('\nWARNING: Some uniprot locations were missing. '
                  + 'Execution will continue with '
                  'default values (missing locations treated as '
                  + self.default_location + '). Make sure to update '
                  'helper file for next execution.\n')

    def cytoplasm_id(self):
        user_id = self._location_map['Cytoplasm']
        if user_id == self._missing_string:
            return 'Cytoplasm'
        else:
            return user_id

    def secreted_id(self):
        user_id = self._location_map['Secreted']
        if user_id == self._missing_string:
            return 'Secreted'
        else:
            return user_id

    def compartment_ids(self):
        ids = []
        for uniprot_id, user_id in self._location_map.iteritems():
            if user_id != self._missing_string:
                ids.append(user_id)
            else:
                ids.append(uniprot_id.replace(' ', '_'))
        return ids
            
    def _read_mapping_data(self):
        """
        Read file containing location id mapping.
        """
        try:
            with open(self._mapping_file, 'rU') as input_stream:
                print('Found file with location mapping data. '
                      'This file will be used to map uniprot locations '
                      'to user-defined locations...')
                # skip header
                next(input_stream)
                # read lines
                for line in input_stream:
                    line = line.rstrip('\n')
                    [location, user_location] = line.split('\t')
                    # check for missing values
                    if (user_location == self._missing_tag) \
                       or (user_location == ''):
                        self._missing_mapping = True
                        user_location = self._missing_string
                    self._location_map[location] = user_location
        except IOError:
            print 'Could not find file with location mapping...'
                
        # update default location (if possible)
        user_default = self._location_map[self.default_location]
        if user_default != self._missing_string:
            self.default_location = user_default

    def _read_curated_data(self):
        """
        Read files containing hand-curated information.
        """
        try:
            with open(self._curation_file, 'rU') as input_stream:
                print('Found file with location data. This file will be used '
                      'to solve ambiguous uniprot annotation...')
                # skip header
                next(input_stream)
                # read lines
                for line in input_stream:
                    line = line.rstrip('\n')
                    [entry, gene_names, name, user_location] = line.split('\t')
                    # check for missing values
                    if (user_location == self._missing_tag) \
                       or (user_location == ''):
                        self._missing_information = True
                    else:
                        # check that user_location is consistent
                        # with known locations
                        if not user_location in self._location_map.values():
                            print('Warning: unknown location ' + user_location
                                  + ' will be replaced by '
                                  + self.default_location
                                  + '. Please check that data given in '
                                  + self._curation_file + ' and '
                                  + self._mapping_file + ' is consistent.')
                        self._curated_locations[entry] \
                            = Location(user_location, gene_names, name)
        except IOError:
            print 'Could not find file with location data...'
            
    def _write_mapping_data(self):
        """
        Write file containing mapping information.
        """
        with open(self._mapping_file, 'w') as output_stream:
            output_stream.write('\t'.join(['UNIPROT NAME', 'USER ID']) + '\n')
            for location, user_location in self._location_map.iteritems():
                if user_location == self._missing_string:
                    user_location = self._missing_tag
                output_stream.write('\t'.join([location, user_location]) + '\n')

    def _write_curated_data(self):
        """
        Write file containing hand-curated information.
        """
        with open(self._curation_file, 'w') as output_stream:
            output_stream.write('\t'.join(['ENTRY', 'GENE NAME', 'NAME',
                                           'LOCATION']) + '\n')
            for entry in self._curated_locations:
                location = self._curated_locations[entry]
                if location.location == self._missing_string:
                    loc = self._missing_tag
                else:
                    loc = location.location
                output_stream.write('\t'.join([entry, location.gene_names,
                                               location.name, loc]) + '\n')
        
    def _extract_uniprot_data(self, uniprot_data):
        """
        Extract location information from uniprot file, using curated data
        if available or writing missing data to file for future curation.
        """
        entry_data = uniprot_data['Entry']
        gene_data = uniprot_data['Gene names']
        name_data = uniprot_data['Protein names']
        location_data = uniprot_data['Subcellular location [CC]']

        location_reader = re.compile(r'SUBCELLULAR LOCATION:\s([\w\s]+\w)')

        curated_data_added = False
        mapping_data_added = False
        for entry, gene_names, name, location_note \
            in zip(entry_data, gene_data, name_data, location_data):
            # if entry is in curated file simply add it
            if self._curated_locations.has_key(entry):
                location = self._curated_locations[entry].location
                # check that the 'curated' data is valid
                if location in self._location_map.values():
                    self.locations[entry] = location
                else:
                    self.locations[entry] = self.default_location
                continue

            # if field is empty, add to data that must be curated
            if pandas.isnull(location_note):
                curated_data_added = True
                self._missing_information = True
                if pandas.isnull(gene_names): gene_names = ''
                if pandas.isnull(name): name = ''
                self._curated_locations[entry] \
                    = Location(self._missing_string, gene_names, name)
                # put default location
                self.locations[entry] = self.default_location
                continue

            # else try to parse field
            try:
                location_name = location_reader.match(location_note).group(1)
            except AttributeError:
                print location_note
                raise
            
            # if location is not in mapping file yet, add it
            if not(self._location_map.has_key(location_name)):
                mapping_data_added = True
                self._missing_mapping = True
                self._location_map[location_name] = self._missing_string

            # return user-defined location
            location = self._location_map[location_name]
            if location != self._missing_string:
                self.locations[entry] = location
            else:
                # if mapping is empty, use uniprot name as compartment name
                self.locations[entry] = location_name.replace(' ', '_')

        if mapping_data_added:
            self._write_mapping_data()
                
        if curated_data_added:
            self._write_curated_data()
