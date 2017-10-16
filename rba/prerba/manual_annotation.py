"""Interface to manual annotation."""

# python 2/3 compatibility
from __future__ import division, print_function, absolute_import

import os.path
from collections import namedtuple

from rba.prerba.curation_data import CurationData


class CurationFile(object):
    """
    Wrapper for curation files.

    Attributes
    ----------
    filename : str
        Location of curation file.
    data : rba.prerba.CurationData
        Curation data.

    """

    def __init__(self, filename, columns):
        """
        Build object from existing curation file or create one.

        Parameters
        ----------
        filename : str
            Location of curation file.
        columns : list of str
            Names of columns of curation file. Names will be overridden if
            a curation file already exists.

        """
        self.filename = filename
        self.data = CurationData(columns)
        try:
            self.data.read(self.filename)
        except IOError:
            print('Helper file {} not found.'.format(filename))


class ManualAnnotation(object):
    """
    Interface to manual annotaton.

    Attributes
    ----------
    subunits : rba.prerba.CurationFile
        Curated protein stoichiometries.
    cofactors : rba.prerba.CurationFile
        Curated protein cofactors.
    locations : rba.prerba.CurationFile
        Curated protein locations.
    location_map : rba.prerba.CurationFile
        Curated map from uniprot locations to user locations.
    unknown_proteins : rba.prerba.CurationFile
        Curated map from sbml genes to uniprot genes.
    metabolites : rba.prerba.CurationFile
        Curated metabolite information.
    macrocomponents : rba.prerba.CurationFile
        Curated macrocomponent information.

    """

    def __init__(self, input_dir):
        """
        Build from input director.

        Parameters
        ----------
        input_dir : str
            Path to directory containing curation files.

        """
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

    def print_warnings(self):
        """Print warning for curation files missing information."""
        if self.subunits.data.has_missing_information():
            print('WARNING: Several uniprot subunit notes were ambiguous. '
                  'Please read file {}. Execution will continue with '
                  'default values.'
                  .format(self.subunits.filename))
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
                  'modify {}.'
                  .format(self.location_map.filename))
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
