"""Find protein information from uniprot and manual data."""

# python 2/3 compatibility
from __future__ import division, print_function, absolute_import

from collections import namedtuple
import pandas

from rba.prerba.uniprot_data import Cofactor, UniprotData
from rba.prerba.manual_annotation import (
    CuratedCofactors, CuratedSubunits, CuratedLocations,
    CuratedLocationMap, CuratedUnknownProteins
    )


class Protein(object):
    """
    Basic protein information.

    Attributes
    ----------
    location : str
        Location of protein.
    stoichiometry : int or float
        Stoichiometry of protein within enzymatic complex.
    cofactors : Cofactor
        Cofactors of protein.
    sequence : str
        Sequence of protein (amino acids in one-lette format).

    """

    def __init__(self):
        """Build default protein."""
        self.location = None
        self.stoichiometry = None
        self.cofactors = None
        self.sequence = None


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
        Map from uniprot locations to user locations.

    """

    def __init__(self, input_dir):
        """
        Build object from uniprot and manual data.

        Parameters
        ----------
        uniprot_data : rba.prerba.uniprot_data.UniprotData
            Uniprot data.
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
        self._default_location = self._location_map.data.get('Cytoplasm',
                                                             'Cytoplasm')
        self._average_id = self.average_protein_id(self._default_location)

    def protein_and_reference(self, gene_id):
        """
        Retrieve base protein information and protein reference.

        Parameters
        ----------
        gene_id : str
            Gene identifier.

        Returns
        -------
        Protein
            Basic protein information. None if protein is unknown or
            'spontaneous'.
        tuple (str, numeric)
            (name, stoichiometry) that should be used in protein reference
            (e.g. for enzymes or other machineries composed of proteins).
            In general, value is (gene_id, protein.stoichiometry).
            If protein is unknown, value is (average_protein, 1).
            If protein is 'spontaneous', None is returned.

        """
        if not gene_id:
            return None, None
        user_id = self._user_ids.data.get(gene_id, 0)
        if user_id != 0:
            if user_id == '' or pandas.isnull(user_id):
                return None, None
            elif user_id.startswith('average_protein_'):
                return None, (user_id, 1)
        protein = self.find_uniprot(user_id if user_id != 0 else gene_id)
        if protein:
            return protein, (gene_id, protein.stoichiometry)
        else:
            if user_id == 0:
                self._user_ids.append(gene_id, self._average_id)
            return None, (self._average_id, 1)

    def find_uniprot(self, gene_id):
        """
        Retrieve information for protein from gene identifier.

        Parameters
        ----------
        gene_id : str
            Gene identifier.

        Returns
        -------
        Protein
            Basic protein information (missing information replaced by
            default values).

        """
        uniprot_id = self._uniprot.entry(gene_id)
        if not uniprot_id:
            return None
        uniprot_line = self._uniprot.line(uniprot_id)
        protein = self._manual_info(uniprot_id)
        if protein.location is None:
            protein.location = self._uniprot_location(uniprot_line)
        if protein.cofactors is None:
            protein.cofactors = self._uniprot_cofactors(uniprot_line)
        if protein.stoichiometry is None:
            protein.stoichiometry = self._uniprot_subunits(uniprot_line)
        protein.sequence = uniprot_line['Sequence']
        return protein

    def average_protein_id(self, compartment):
        """Return identifier of average protein in given compartment."""
        return 'average_protein_' + compartment

    def average_composition(self):
        """Return average composition of proteins."""
        return self._uniprot.average_protein_composition()

    def compartments(self):
        """
        Return list of compartment identifiers.

        Returns
        -------
        list
            List of compartment identifiers.

        """
        return list(set(self._location_map.data.values()))

    def compartment(self, compartment):
        """
        Return user identifier associated with compartment..

        Parameters
        ----------
        compartment: str
            Valid uniprot compartment identifier (e.g. Cytoplasm, Secreted)
            or user identifier.

        Returns
        -------
        str
            User identifier associated with compartment. If a user identifier
            was provided as an input, it is returned without change.

        """
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

    def _manual_info(self, identifier):
        protein = Protein()
        protein.location = self._locations.data.get(identifier)
        protein.cofactors = self._cofactors.data.get(identifier)
        protein.stoichiometry = self._subunits.data.get(identifier)
        return protein

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
            self._subunits.append(uniprot_line, 1)
            subunits = 1
        return subunits

    def _check_user_identifiers(self):
        invalid_identifiers = [i for i in self._user_ids.data.values()
                               if i and not i.startswith('average_protein_')
                               and not self._uniprot.entry(i)]
        if invalid_identifiers:
            print('Warning: {} are invalid gene identifiers. '
                  'Check data provided provided in {}.'
                  .format(', '.join(invalid_identifiers),
                          self._user_ids._raw_data.filename))

    def _check_location_validity(self):
        known_locations = set(self._location_map.data.keys()
                              + self._location_map.data.values())
        invalid_locations = []
        for loc in set(self._locations.data.values()):
            if not pandas.isnull(loc) and loc not in known_locations:
                invalid_locations.append(loc)
        if invalid_locations:
            print('Warning: unknown location(s) {} will be replaced by {}.'
                  'Check that data in {} and {} are consistent.'
                  .format(', '.join(set(invalid_locations)),
                          self.default_location,
                          self._manual.location_map.filename,
                          self._manual.locations.filename))
