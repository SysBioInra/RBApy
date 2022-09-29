"""Interface to manual annotation."""

# python 2/3 compatibility
from __future__ import division, print_function, absolute_import

import os.path
import warnings
from collections import namedtuple
import pandas

from rba.prerba.curation_data import CurationData
from rba.prerba.uniprot_data import Cofactor

Metabolite = namedtuple('Metabolite', 'name sbml_id concentration')


class CuratedData(object):
    def __init__(self, filename, columns):
        self._raw_data = CurationData(filename, columns)
        self.data = {}
        self._warning = ''

    def update_file(self):
        if self._raw_data.update_file() and self._warning:
            warnings.warn(self._warning, UserWarning)


class CuratedSubunits(CuratedData):
    def __init__(self, input_dir):
        filename = os.path.join(input_dir, 'subunits.tsv')
        super(CuratedSubunits, self).__init__(
            filename,
            ['ENTRY', 'STOICHIOMETRY', 'GENE NAMES',
             'PROTEIN NAME', 'UNIPROT NOTE']
            )
        if self._raw_data.has_missing_information('STOICHIOMETRY'):
            raise UserWarning(filename + ': please fill in the'
                              ' STOICHIOMETRY column.')
        self.data = {r[0]: float(r[1]) for r in self._raw_data.rows()}
        self._warning = (
            'WARNING: ambiguous subunit notes have been added to '
            'file {}. Execution will continue with default values.'
            .format(filename)
            )

    def append(self, uniprot_line, value):
        """Add ambiguous stoichiometry."""
        self.data[uniprot_line.name] = value
        self._raw_data.add_row(
            (uniprot_line.name, value) +
            tuple(uniprot_line[['Gene Names', 'Protein names',
                                'Subunit structure']])
            )


class CuratedLocations(CuratedData):
    def __init__(self, input_dir):
        filename = os.path.join(input_dir, 'locations.tsv')
        super(CuratedLocations, self).__init__(
            filename, ['ENTRY', 'GENE NAME', 'PROTEIN NAME', 'LOCATION']
            )
        self.data = {r[0]: r[3] for r in self._raw_data.rows()}
        if self._raw_data.has_missing_information('LOCATION'):
            raise UserWarning(filename + ': please fill in the'
                              ' LOCATION column.')
        self._warning = (
            'WARNING: ambiguous UniProt locations have been added to '
            'file {}. Execution will continue with default values.'
            .format(filename)
            )

    def append(self, uniprot_line, value):
        """Add ambiguous location."""
        self.data[uniprot_line.name] = value
        self._raw_data.add_row(
            (uniprot_line.name,) +
            tuple(uniprot_line[['Gene Names', 'Protein names']]) +
            (value,)
            )


class CuratedCofactors(CuratedData):
    def __init__(self, input_dir):
        filename = os.path.join(input_dir, 'cofactors.tsv')
        super(CuratedCofactors, self).__init__(
            filename,
            ['ENTRY', 'CHEBI', 'NAME',
             'STOICHIOMETRY', 'UNIPROT ANNOTATION']
            )
        self.data = {}
        for row in self._raw_data.rows():
            self._add_to_data(row[0], Cofactor(*row[1:]))
        self._warning = (
            'WARNING: ambiguous UniProt cofactor notes have been added to '
            'file {}. Execution will continue with default values.'
            .format(filename, CurationData.missing_tag)
            )

    def append(self, entry, cofactors):
        """Add ambiguous cofactors."""
        for c in cofactors:
            self._add_to_data(entry, c)
        self._raw_data.add_rows([(entry,) + c for c in cofactors])

    def _add_to_data(self, entry, cofactor):
        """Add cofactor to data only if it has valid information."""
        # do not include cofactors with missing chebi or stoichiometry
        # as well as cofactors with 0 stoichiometry
        sto = cofactor.stoichiometry
        cofactor_list = self.data.setdefault(entry, [])
        if (pandas.notnull(cofactor.chebi)
                and pandas.notnull(sto)
                and float(sto) > 0):
            cofactor_list.append(cofactor)


class CuratedLocationMap(CuratedData):
    def __init__(self, input_dir):
        filename = os.path.join(input_dir, 'location_map.tsv')
        super(CuratedLocationMap, self).__init__(
            filename, ['UNIPROT NAME', 'USER ID']
            )
        if self._raw_data.has_missing_information('USER ID'):
            raise UserWarning(filename + ': please fill in the'
                              ' USER ID column.')
        self.data = {r[0]: r[1] for r in self._raw_data.rows()}
        # add mandatory compartments (if they are missing)
        self.data.setdefault('Secreted', 'Secreted')
        self._warning = (
            'WARNING: UniProt locations with no user-defined '
            'counterpart have been added to {}.'
            .format(filename)
            )

    def append(self, location, default_value):
        """Add location without user counterpart."""
        self.data[location] = default_value
        self._raw_data.add_row((location, default_value))


class CuratedUnknownProteins(CuratedData):
    def __init__(self, input_dir):
        filename = os.path.join(input_dir, 'unknown_proteins.tsv')
        super(CuratedUnknownProteins, self).__init__(
            filename, ['SBML ID', 'UNIPROT GENE']
            )
        if self._raw_data.has_missing_information('UNIPROT GENE'):
            raise UserWarning(filename + ': please fill in the'
                              ' UNIPROT GENE column.')
        self.data = {r[0]: r[1] for r in self._raw_data.rows()}
        self._warning = (
            'WARNING: SBML genes not referenced in UniProt have been added to '
            'file {}. Execution will continue with default values.'
            .format(filename)
            )

    def append(self, gene_id, default_value):
        self.data[gene_id] = default_value
        self._raw_data.add_row((gene_id, default_value))


class CuratedMetabolites(CuratedData):
    def __init__(self, input_dir, known_species):
        filename = os.path.join(input_dir, 'metabolites.tsv')
        super(CuratedMetabolites, self).__init__(
            filename, ['ID', 'NAME', 'SBML ID', 'CONCENTRATION']
            )
        self.data = {}
        invalid_ids = []
        for id_, name, sbml_id, conc in self._raw_data.rows():
            if pandas.isnull(sbml_id):
                sbml_id = None
            elif sbml_id not in known_species:
                invalid_ids.append(id_)
                sbml_id = None
            if pandas.isnull(conc) or conc == '':
                conc = 0
            self.data[id_] = Metabolite(name, sbml_id, float(conc))
        if invalid_ids:
            print(
                '{}: {} do not map to valid SBML metabolite ids. '
                'These metabolites will be removed from processes and '
                'production targets.'
                .format(filename, ', '.join(invalid_ids))
                )
        self._warning = (
            'WARNING: unidentified key metabolites have been added to file {}.'
            ' These metabolites will be removed from processes and '
            'production targets.'
            .format(filename)
            )

    def append(self, id_, name, sbml_id, concentration):
        """Add unrecognized metabolite."""
        self.data[id_] = Metabolite(name, sbml_id, float(concentration))
        self._raw_data.add_row((id_, name, sbml_id, concentration))


class CuratedMacrocomponents(CuratedData):
    def __init__(self, input_dir, known_species):
        filename = os.path.join(input_dir, 'macrocomponents.tsv')
        super(CuratedMacrocomponents, self).__init__(
            filename, ['TARGET_METABOLITE', 'TARGET_FLUX']
            )
        self.data = {}
        invalid_ids = []
        for met, flux in self._raw_data.rows():
            if met in known_species:
                self.data[met] = float(flux)
            else:
                invalid_ids.append(met)
        if invalid_ids:
            print(
                '{}: {} are not valid SBML metabolite ids. '
                'These entries will be removed from production targets.'
                .format(filename, ', '.join(invalid_ids))
                )
