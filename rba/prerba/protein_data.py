"""Find protein information from UniProt and manual data."""

# python 2/3 compatibility
from __future__ import division, print_function, absolute_import

from itertools import chain
import pandas

from rba.prerba.macromolecule import Protein
from rba.prerba.uniprot_data import UniprotData
from rba.prerba.manual_annotation import (
    CuratedCofactors, CuratedSubunits, CuratedLocations,
    CuratedLocationMap, CuratedUnknownProteins
    )


class ProteinData(object):
    """
    Interface to protein data.

    Attributes
    ----------
    default_stoichiometry : int or float
        Default protein stoichiometry.
    default_location : str
        Default protein location.
    location_map : dict
        Map from UniProt locations to user locations.

    """
    def __init__(self, input_dir):
        """
        Build object from UniProt and manual data.

        Parameters
        ----------
        uniprot_data : rba.prerba.uniprot_data.UniprotData
            UniProt data.
        manual_data : rba.prerba.manual_annotation.ManualAnnotation
            Curated data.

        """
        self._uniprot = UniprotData(input_dir)
        self._locations = CuratedLocations(input_dir)
        self._cofactors = CuratedCofactors(input_dir)
        self._subunits = CuratedSubunits(input_dir)
        self._user_ids = CuratedUnknownProteins(input_dir)
        self._location_map = CuratedLocationMap(input_dir)
        self._check_location_validity()
        self._check_user_identifiers()
        self._default_location = self._location_map.data.get(
            'Cytoplasm', 'Cytoplasm'
            )
        self._average_id = self.average_protein_id(self._default_location)

    def _check_location_validity(self):
        known_locations = set(chain(self._location_map.data.keys(),
                                    self._location_map.data.values()))
        invalid_locations = []
        for loc in set(self._locations.data.values()):
            if not pandas.isnull(loc) and loc not in known_locations:
                invalid_locations.append(loc)
        self._warn_invalid_locations(invalid_locations)

    def _warn_invalid_locations(self, invalid_locations):
        if invalid_locations:
            print('Warning: unknown location(s) {} will be replaced by {}.'
                  'Check that data in {} and {} are consistent.'
                  .format(', '.join(set(invalid_locations)),
                          self._default_location,
                          self._manual.location_map.filename,
                          self._manual.locations.filename))

    def _check_user_identifiers(self):
        invalid_identifiers = [i for i in self._user_ids.data.values()
                               if i and not i.startswith('average_protein_')
                               and not self._uniprot.entry(i)]
        if invalid_identifiers:
            print('Warning: {} are invalid gene identifiers. '
                  'Check data provided provided in {}.'
                  .format(', '.join(invalid_identifiers),
                          self._user_ids._raw_data.filename))

    def create_protein_from_gene_id(self, gene_id):
        """
        Retrieve base protein information.

        Returns
        -------
        Protein
            Basic protein information. None if protein is unknown,
            average protein or 'spontaneous'.
        """
        user_id = self._user_ids.data.get(gene_id, gene_id)
        if (self._is_spontaneous_id(user_id)
                or self._is_average_protein_id(user_id)):
            return None
        uniprot_id = self._uniprot.entry(user_id)
        if uniprot_id:
            protein = self.create_protein_from_uniprot_id(uniprot_id)
            protein.id = gene_id
            return protein
        else:
            # unknown gene identifier
            self._user_ids.append(gene_id, self._average_id)
            return None

    def _is_spontaneous_id(self, id_):
        return id_ == '' or pandas.isnull(id_)

    def _is_average_protein_id(self, id_):
        return(id_.startswith('average_protein_') and
               id_[16:] in self._location_map.data.values())

    def create_protein_from_uniprot_id(self, uniprot_id):
        protein = Protein()
        self._fill_with_manual_info(protein, uniprot_id)
        self._fill_with_uniprot_info(protein, uniprot_id)
        return protein

    def _fill_with_manual_info(self, protein, uniprot_id):
        protein.location = self._locations.data.get(uniprot_id)
        protein.cofactors = self._cofactors.data.get(uniprot_id)
        protein.stoichiometry = self._subunits.data.get(uniprot_id)

    def _fill_with_uniprot_info(self, protein, uniprot_id):
        uniprot_line = self._uniprot.line(uniprot_id)
        if protein.location is None:
            protein.location = self._uniprot_location(uniprot_line)
        if protein.cofactors is None:
            protein.cofactors = self._uniprot_cofactors(uniprot_line)
        if protein.stoichiometry is None:
            protein.stoichiometry = self._uniprot_subunits(uniprot_line)
        if not protein.sequence:
            protein.sequence = uniprot_line['Sequence']

    def _uniprot_location(self, uniprot_line):
        location = self._uniprot.find_location(uniprot_line)
        if location:
            user_location = self._location_map.data.get(location)
            if not user_location:
                default_name = location.replace(' ', '_')
                self._location_map.append(location, default_name)
                location = default_name
            else:
                location = user_location
        else:
            self._locations.append(uniprot_line, self._default_location)
            location = self._default_location
        return location

    def _uniprot_cofactors(self, uniprot_line):
        cofactors, curation_needed \
            = self._uniprot.find_cofactors(uniprot_line)
        if curation_needed:
            self._cofactors.append(uniprot_line.name, cofactors)
            cofactors = self._cofactors.data.get(uniprot_line.name, [])
        return cofactors

    def _uniprot_subunits(self, uniprot_line):
        subunits = self._uniprot.find_subunits(uniprot_line)
        if not subunits:
            subunits = 1
        self._subunits.append(uniprot_line, subunits)
        return subunits

    def reference(self, gene_id):
        """
        Retrieve protein reference.

        Returns
        -------
        tuple (str, numeric)
            (id, stoichiometry) that should be used in protein reference.
            If protein is 'spontaneous', None is returned.

        """
        user_id = self._user_ids.data.get(gene_id, gene_id)
        if self._is_spontaneous_id(user_id):
            return None
        if self._is_average_protein_id(user_id):
            return (user_id, 1)
        uniprot_id = self._uniprot.entry(gene_id)
        if uniprot_id:
            stoichiometry = self._stoichiometry(uniprot_id)
            return (gene_id, stoichiometry)
        else:
            # unknown gene identifier
            self._user_ids.append(gene_id, self._average_id)
            return (self._average_id, 1)

    def _stoichiometry(self, uniprot_id):
        result = self._subunits.data.get(uniprot_id)
        if result:
            return result
        else:
            return self._uniprot_subunits(self._uniprot.line(uniprot_id))

    def average_protein_id(self, compartment):
        """Return identifier of average protein in given compartment."""
        return 'average_protein_' + compartment

    def average_composition(self):
        """Return average composition of proteins."""
        return self._uniprot.average_protein_composition()

    def compartments(self):
        """Return list of compartment identifiers."""
        return list(set(self._location_map.data.values()))

    def compartment(self, compartment):
        """Return user identifier associated with compartment."""
        if compartment in self._location_map.data.values():
            return compartment
        else:
            return self._location_map.data[compartment]

    def update_helper_files(self):
        """Update helper files (if needed)."""
        self._locations.update_file()
        self._location_map.update_file()
        self._cofactors.update_file()
        self._subunits.update_file()
        self._user_ids.update_file()
