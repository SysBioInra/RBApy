"""Interface to manual data."""

# python 2/3 compatibility
from __future__ import division, print_function

from collections import namedtuple
import pandas

from rba.prerba.uniprot_data import Cofactor


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
        self.stoichiometry = 0
        self.cofactors = []
        self.sequence = None


class AmbiguousData(object):
    """
    Ambiguous protein information.

    Attributes
    ----------
    locations : list of tuples
        Ambiguous locations in curation file format.
    cofactors : list of tuples
        Ambiguous cofactors in curation file format.
    unknow_locations : list of tuples
        Location without user counterpart in curation file format.
    subunits : list of tuples
        Ambiguous stoichiometries in curation file format.

    """

    def __init__(self):
        """Build default object."""
        self.reset()

    def reset(self):
        """Remove all currently stored information."""
        self.locations = []
        self.cofactors = []
        self.unknown_locations = []
        self.subunits = []

    def add_location(self, uniprot_line, default_value):
        """
        Add ambiguous location.

        Parameters
        ----------
        uniprot_line : pandas.Series
            Raw uniprot data.
        default_value : str
            Default location.

        """
        self.locations.append(
            (uniprot_line.name,) +
            tuple(uniprot_line[['Gene names', 'Protein names']]) +
            (default_value,)
            )

    def add_unknown_location(self, location, default_value):
        """
        Add location without user counterpart.

        Parameters
        ----------
        location : str
            Unknown location.
        default_value : str
            Default location.

        """
        self.unknown_locations.append((location, default_value))

    def add_cofactors(self, entry, cofactors):
        """
        Add ambiguous cofactors.

        Parameters
        ----------
        entry : str
            Uniprot identifier of protein.
        cofactors : list of Cofactor
            Cofactors containing ambiguous information.

        """
        self.cofactors += [(entry,) + c for c in cofactors]

    def add_subunits(self, uniprot_line, default_value):
        """
        Add ambiguous stoichiometry.

        Parameters
        ----------
        uniprot_line : pandas.Series
            Raw uniprot data.
        default_value : str
            Default stoichiometry.

        """
        self.subunits.append(
            (uniprot_line.name, default_value) +
            tuple(uniprot_line[['Gene names', 'Protein names',
                                'Subunit structure [CC]']])
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
        Map from uniprot locations to user locations.

    """

    def __init__(self, uniprot_data, manual_data):
        """
        Build object from uniprot and manual data.

        Parameters
        ----------
        uniprot_data : rba.prerba.uniprot_data.UniprotData
            Uniprot data.
        manual_data : rba.prerba.manual_annotation.ManualAnnotation
            Curated data.

        """
        self.default_stoichiometry = 1
        self.default_location = 'Cytoplasm'
        self._uniprot = uniprot_data
        # import manual data
        self._subunits = {
            row[0]: row[1] for row in manual_data.subunits.data.data.values
            }
        self._cofactors = {}
        for row in manual_data.cofactors.data.data.values:
            self._cofactors.setdefault(row[0], []).append(Cofactor(*row[1:]))
        self._locations = {
            row[0]: row[3] for row in manual_data.locations.data.data.values
        }
        self.location_map = {
            row[0]: row[1] for row in manual_data.location_map.data.data.values
            }
        # update default location (if possible)
        user_default = self.location_map.get(self.default_location, None)
        if not pandas.isnull(user_default):
            self._default_location = user_default
        # check that all locations are valid
        known_locations = set(self.location_map.values())
        for loc in set(self._locations.values()):
            if not pandas.isnull(loc) and loc not in known_locations:
                print('Warning: unknown location {} will be replaced by {}.'
                      'Check that data in {} and {} is consistent.'
                      .format(location, self.default_location,
                              manual_data.location_map.filename,
                              manual_data.locations.filename))
        self._ambiguous = AmbiguousData()

    def find_protein(self, uniprot_id):
        """
        Retrieve information for protein with given id.

        Parameters
        ----------
        uniprot_id : str
            Uniprot identifier.

        Returns
        -------
        Protein
            Basic protein information (missing information replaced by
            default values).

        """
        NOT_FOUND = -1
        protein = Protein()
        uniprot_line = self._uniprot.line(uniprot_id)
        # location
        location = self._locations.get(uniprot_id, NOT_FOUND)
        if location == NOT_FOUND:
            location = self._uniprot.find_location(uniprot_line)
            if location:
                user_location = self.location_map.get(location, None)
                if user_location:
                    location = user_location
                else:
                    standard_loc = location.replace(' ', '_')
                    self.location_map[location] = standard_loc
                    self._ambiguous.add_unknown_location(location,
                                                         standard_loc)
                    location = standard_loc
            else:
                location = self.default_location
                self._locations[uniprot_id] = None
                self._ambiguous.add_location(uniprot_line, None)
        protein.location = location
        # cofactors
        cofactors = self._cofactors.get(uniprot_id, None)
        if not cofactors:
            cofactors, to_cure = self._uniprot.find_cofactors(uniprot_line)
            self._cofactors[uniprot_id] = cofactors
            self._ambiguous.add_cofactors(uniprot_id, to_cure)
        protein.cofactors = cofactors
        # subunits
        subunits = self._subunits.get(uniprot_id, None)
        if not subunits:
            subunits = self._uniprot.find_subunits(uniprot_line)
            if not subunits:
                subunits = self.default_stoichiometry
                self._subunits[uniprot_id] = self.default_stoichiometry
                self._ambiguous.add_subunits(uniprot_line,
                                             self.default_stoichiometry)
        protein.stoichiometry = subunits
        protein.sequence = uniprot_line['Sequence']
        self.fill_missing_values(protein)
        return protein

    def fill_missing_values(self, protein):
        """
        Replace missing information with default values.

        Parameters
        ----------
        protein : Protein
            Protein information. Missing values will be replaced by default
            values in place.

        """
        if pandas.isnull(protein.location):
            protein.location = self.default_location
        if pandas.isnull(protein.stoichiometry):
            protein.stoichiometry = self.default_stoichiometry
        # remove cofactors with unknown chebi or zero stoichiometry
        # replace missing stoichiometry with default stoichiometry
        cofactors = []
        for cofactor in protein.cofactors:
            if pandas.notnull(cofactor.chebi) and cofactor.stoichiometry != 0:
                if pandas.isnull(cofactor.stoichiometry):
                    cofactors.append(cofactor._replace(
                        stoichiometry=self.default_stoichiometry
                        ))
                else:
                    cofactors.append(cofactor)
        protein.cofactors = cofactors

    def update_manual_annotation(self, manual_data):
        """
        Add all ambiguous data previously retrieved to curation files.

        Parameters
        ----------
        manual_data : rba.prerba.manual_annotation.ManualAnnotation
            Curation files where data should be added to.

        """
        update_and_write(manual_data.locations, self._ambiguous.locations)
        update_and_write(manual_data.subunits, self._ambiguous.subunits)
        update_and_write(manual_data.cofactors, self._ambiguous.cofactors)
        update_and_write(manual_data.location_map,
                         self._ambiguous.unknown_locations)
        self._ambiguous.reset()
        manual_data.print_warnings()


def update_and_write(curation_file, new_data):
    """Write data to curation file if necessary."""
    if new_data:
        curation_file.data.add(new_data)
        curation_file.data.write(curation_file.filename)
