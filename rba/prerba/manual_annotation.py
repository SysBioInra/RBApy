"""Interface to manual annotation."""

# python 2/3 compatibility
from __future__ import division, print_function, absolute_import

import os.path
from collections import namedtuple

from rba.prerba.curation_data import CurationData

class CurationFile(object):
    def __init__(self, filename, columns):
        self.filename = filename
        self.data = CurationData(columns)
        try:
            self.data.read(self.filename)
        except IOError:
            print('Helper file {} not found.'.format(filename))


class ManualAnnotation(object):
    """
    Interface to manual annotaton.
    """

    def __init__(self, input_dir):
        self.subunits = CurationFile(
            os.path.join(input_dir, 'subunits.tsv'),
            ['ENTRY', 'STOICHIOMETRY', 'GENE NAMES',
             'PROTEIN NAME', 'UNIPROT NOTE']
            )
        self.cofactors = CurationFile(
            os.path.join(input_dir, 'cofactors.tsv'),
            ['ENTRY', 'CHEBI', 'NAME',
             'STOICHIOMETRY', 'UNIPROT ANNOTATION']
            )
        self.locations = CurationFile(
            os.path.join(input_dir, 'locations.tsv'),
            ['ENTRY', 'GENE NAME', 'NAME', 'LOCATION']
            )
        self.location_map = CurationFile(
            os.path.join(input_dir, 'location_map.tsv'),
            ['UNIPROT NAME', 'USER ID']
            )
        self.unknown_proteins = CurationFile(
            os.path.join(input_dir, 'unknown_proteins.tsv'),
            ['SBML ID', 'UNIPROT GENE']
            )
        self.metabolites = CurationFile(
            os.path.join(input_dir, 'metabolites.tsv'),
            ['ID', 'NAME', 'SBML ID', 'CONCENTRATION']
            )
        self.macrocomponents = CurationFile(
            os.path.join(input_dir, 'macrocomponents.tsv'),
            ['TARGET_METABOLITE', 'TARGET_FLUX']
            )

    def write_files(self):
        for curation in [self.subunits, self.location_map,
                         self.cofactors, self.locations]:
            curation.data.write(curation.filename)

    def print_warnings(self):
        if self.subunits.data.has_missing_information():
            print('WARNING: Several uniprot notes were ambiguous. '
                  'Please read file {}, check data and specify all '
                  'information tagged as {}. Execution will continue with '
                  'default values.'
                  .format(self.subunits.filename,
                          CurationData.missing_tag))
        if self.cofactors.data.has_missing_information():
            print('WARNING: Several uniprot notes were ambiguous. '
                  'Please read file {}, check data and specify all '
                  'information tagged as {}. Add lines if necessary but avoid '
                  'removing lines. Execution will continue with default '
                  'values.'
                  .format(self.cofactors.filename,
                          CurationData.missing_tag))
        if self.location_map.data.has_missing_information():
            print('WARNING: Some uniprot locations had no user-defined '
                  'counterpart. Uniprot names will be used. '
                  'If you wish to change the compartment names, please '
                  'read file {} and specify all information tagged as {}.'
                  .format(self.location_map.filename,
                          CurationData.missing_tag))
        if self.locations.data.has_missing_information():
            print('WARNING: Some uniprot locations were missing. '
                  'Execution will continue with default value. '
                  'You can specify missing locations in file {}.'
                  .format(self.locations.filename))
        if self.unknown_proteins.data.has_missing_information():
            print('WARNING: Several enzymatic proteins defined in your SBML '
                  'could not be retrieved. '
                  'Please read file {}, check data and specify all missing '
                  'information. Execution will continue with default '
                  'values (unknown proteins replaced by average enzymatic '
                  'protein in the cytosol).'
                  .format(self.unknown_proteins.filename))
        if self.metabolites.data.has_missing_information():
            print('WARNING: Several key metabolites could not be identified. '
                  'Please read file {}, check data and specify all missing '
                  'information. Unknown metabolites will be removed '
                  'from production targets for this run.'
                  .format(self.metabolites.filename))
