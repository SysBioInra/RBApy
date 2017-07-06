"""
Module defining LocationData and LocationParser classes.
"""

# python 2/3 compatibility imports
from __future__ import division, print_function

# global imports
import os.path
import re
import pandas

# local imports
from rba.prerba.curation_data import CurationData

class LocationData(object):
    """
    Class used to parse uniprot location field.

    Attributes:
        locations: dict mapping uniprot entries with their compartment.
    """

    _mapping_file = 'location_map.tsv'
    _helper_file = 'locations.tsv'
    def __init__(self, uniprot_data, input_dir='.'):
        """
        Constructor from pandas.DataFrame.

        Args:
            uniprot_data: uniprot data as a pandas.DataFrame.
            input_dir: path to folder containing helper files.
        """
        self.default_location = 'Cytoplasm'

        print('\nParsing locations'
              '\n-----------------')

        # load mapping data
        mapping_file = os.path.join(input_dir, self._mapping_file)
        mapping_data = CurationData(['UNIPROT NAME', 'USER ID'])
        try:
            mapping_data.read(mapping_file)
            print('Mapping file found.')
        except IOError:
            print('Mapping file not found.')
        self._location_map = {x[0]: x[1] for x in mapping_data.data.values}
        # update default location (if possible)
        user_default = self._location_map.get(self.default_location, None)
        if not pandas.isnull(user_default):
            self.default_location = user_default

        # load curated data
        curation_file = os.path.join(input_dir, self._helper_file)
        curated_data = CurationData(['ENTRY', 'GENE NAME', 'NAME', 'LOCATION'])
        try:
            curated_data.read(curation_file)
            print('Helper file found.')
        except IOError:
            print('Helper file not found.')
        curated_locations = {}
        known_locations = list(self._location_map.values())
        for row in curated_data.data.values:
            location = row[3]
            curated_locations[row[0]] = location
            if not pandas.isnull(location) and location not in known_locations:
                print('Warning: unknown location {} will be replaced by {}.'
                      'Please check that data given in {} and {} is consistent.'
                      .format(location, self.default_location,
                              curation_file, mapping_file))

        # read uniprot data
        self.locations, data_to_cure, mapping_to_cure \
            = self._extract_data(uniprot_data, curated_locations,
                                 self._location_map)
        if mapping_to_cure:
            mapping_data.add(mapping_to_cure)
            mapping_data.write(mapping_file)
            self._location_map.update(mapping_to_cure)
        if data_to_cure:
            curated_data.add(data_to_cure)
            curated_data.write(curation_file)

        # print warning if information is missing
        if mapping_data.has_missing_information():
            print('WARNING: Some uniprot locations had no user-defined '
                  'counterpart. Uniprot names will be used.'
                  'If you wish to change the compartment names, please '
                  'read file {} and specify all information tagged as {}.'
                  .format(mapping_file, CurationData.missing_tag))
        if curated_data.has_missing_information():
            print('WARNING: Some uniprot locations were missing. '
                  'Execution will continue with default value ({}). '
                  'You can specify missing locations in file {}.'
                  .format(self.default_location, curation_file))

    def cytoplasm_id(self):
        """
        Accessor to compartment representing cytoplasm.

        Returns:
            Identifier of cytoplasm.
        """
        user_id = self._location_map['Cytoplasm']
        return 'Cytoplasm' if pandas.isnull(user_id) else user_id

    def secreted_id(self):
        """
        Accessor to compartment representing extracellular space.

        Returns:
            Identifier of extracellular space.
        """
        user_id = self._location_map['Secreted']
        return 'Secreted' if pandas.isnull(user_id) else user_id

    def compartment_ids(self):
        """
        Accessor to all compartment identifiers.

        Returns:
            List of compartment identifiers.
        """
        ids = []
        for uniprot_id, user_id in self._location_map.items():
            if not pandas.isnull(user_id):
                ids.append(user_id)
            else:
                ids.append(uniprot_id.replace(' ', '_'))
        return ids

    def _extract_data(self, uniprot_data, curated_locations, location_map):
        """
        Extract location information from uniprot file.
        """
        entry_data = uniprot_data['Entry']
        gene_data = uniprot_data['Gene names']
        name_data = uniprot_data['Protein names']
        location_data = uniprot_data['Subcellular location [CC]']
        parser = LocationParser()

        locations = {}
        data_to_cure = []
        mapping_to_cure = {}
        known_locations = list(location_map.values())
        for entry, gene_names, name, location_note \
            in zip(entry_data, gene_data, name_data, location_data):
            # if entry is in curated file add it, otherwise parse field
            if entry in curated_locations:
                location = curated_locations[entry]
                # check that the 'curated' data is valid
                if location in known_locations:
                    locations[entry] = location
                else:
                    locations[entry] = self.default_location
            elif pandas.isnull(location_note):
                # if field is empty, put default location
                locations[entry] = self.default_location
                data_to_cure.append((entry, gene_names, name, None))
            else:
                # parse field
                location_name = parser.parse(location_note)
                # if location is not in mapping file or mapping was not
                # defined, use uniprot name as compartment name
                if location_name not in location_map:
                    mapping_to_cure[location_name] = None
                    location = None
                else:
                    location = self._location_map[location_name]
                if pandas.isnull(location):
                    locations[entry] = location_name.replace(' ', '_')
                else:
                    locations[entry] = location
        return locations, data_to_cure, list(mapping_to_cure.items())


class LocationParser(object):
    """
    Class parsing 'Subcellular location' field of uniprot.
    """

    _location_reader = re.compile(r'SUBCELLULAR LOCATION:\s([\w\s]+\w)')

    def parse(self, field):
        """
        Parse 'Subcellular location' field in uniprot.

        Returns:
            Compartment read.
        """
        try:
            return self._location_reader.match(field).group(1)
        except AttributeError:
            print(field)
            raise
