"""Module defining UniprotData class."""

# python 2/3 compatibility
from __future__ import division, print_function, absolute_import

# global imports
from collections import Counter, namedtuple
import os.path
import re
import pandas

Cofactor = namedtuple('Cofactor', 'chebi name stoichiometry uniprot_note')


class UniprotData(object):
    """
    Class parsing RBA-relevant Uniprot data.

    Parameters
    ----------
    data: pandas.DataFrame
        raw uniprot data

    """

    def __init__(self, input_dir):
        """
        Build from input directory.

        Parameters
        ----------
        input_file: path to uniprot file.

        """
        # open uniprot data
        self.data = pandas.read_csv(os.path.join(input_dir, 'uniprot.csv'),
                                    sep='\t')
        self.data.set_index('Entry', inplace=True)
        # create mapping from gene ids to uniprot ids
        self._gene_to_entry = {}
        gene_reader = re.compile(r'([^\s]+)')
        for entry, genes in zip(self.data.index, self.data['Gene names']):
            # transform raw uniprot field into standardized list
            if pandas.isnull(genes):
                continue
            gene_ids = set(g.upper() for g in gene_reader.findall(genes))
            for gene in gene_ids:
                self._gene_to_entry[gene] = entry
        # create parsers
        self._location_parser = LocationParser()
        self._cofactor_parser = CofactorParser()
        self._subunit_parser = SubunitParser()

    def line(self, uniprot_id):
        """
        Return data line corresponding to uniprot identifier.

        Parameters
        ----------
        uniprot_id : str
            Uniprot identifier of a protein.

        Returns
        -------
        pandas.Series
            Data associated with protein.

        """
        return self.data.loc[uniprot_id]

    def find_location(self, uniprot_line):
        """
        Parse location of protein.

        Parameters
        ----------
        uniprot_line : pandas.Series
            Protein data.

        Returns
        -------
        str
            Standardized location of protein.

        """
        return self._location_parser.parse(
            uniprot_line['Subcellular location [CC]']
            )

    def find_cofactors(self, uniprot_line):
        """
        Parse cofactors of protein.

        Parameters
        ----------
        uniprot_line : pandas.Series
            Protein data.

        Returns
        -------
        str
            Standardized cofactors of protein.

        """
        return self._cofactor_parser.parse(uniprot_line['Cofactor'])

    def find_subunits(self, uniprot_line):
        """
        Parse stoichiometry of protein.

        Parameters
        ----------
        uniprot_line : pandas.Series
            Protein data.

        Returns
        -------
        str
            Standardized stoichiometry of protein.

        """
        return self._subunit_parser.parse(
            uniprot_line['Subunit structure [CC]']
            )

    def entry(self, gene):
        """
        Find uniprot entries from gene identifiers.

        Parameters
        ----------
        gene_ids : list
            Name of genes to retrieve.

        Returns
        -------
        result : dict
            Dictionary where keys are gene ids and values are
            corresponding uniprot entries.
        not_found : list
            Gene ids that could not be retrieved.

        """
        return self._gene_to_entry.get(gene.upper(), None)

    def average_protein_composition(self):
        """
        Compute average protein composition.

        Returns
        -------
        dict
            Dictionary where keys are amino acids (one letter format) and
            values their average number in a protein.

        """
        composition = Counter()
        for sequence in self.data['Sequence']:
            composition.update(sequence)
        nb_proteins = len(self.data.index)
        for aa in composition:
            composition[aa] /= nb_proteins
        return dict(composition)


class LocationParser(object):
    """Class parsing 'Subcellular location' field of uniprot."""

    _location_reader = re.compile(r'SUBCELLULAR LOCATION:\s([\w\s]+\w)')

    def parse(self, field):
        """
        Parse 'Subcellular location' field in uniprot.

        Parameters
        ----------
        field : str
            Subcellular location field from uniprot.

        Returns
        -------
        str
            Compartment read.

        """
        if pandas.isnull(field):
            return None
        try:
            return self._location_reader.match(field).group(1)
        except AttributeError:
            print(field)
            raise


class SubunitParser(object):
    """
    Class parsing 'Subunit' uniprot field.

    Attributes
    ----------
    prefix_rule : dict
        Dictionary determining rule used to infer stoichiometry.
        Keys are all caps prefixes preceding 'mer' in words found
        in uniprot field,
        values are stoichiometries associated with them. For example,
        prefix_rule[MONO] = 1.

    """

    prefix_rule = {'MONO': 1, 'HETERODI': 1, 'HOMODI': 2, 'HOMOTRI': 3,
                   'HOMOTETRA': 4, 'HOMOPENTA': 5, 'HOMOHEXA': 6,
                   'HEPTA': 7, 'HOMOOCTA': 8, 'HOMODECA': 10, 'HOMODODECA': 12}
    _subunit_reader = re.compile(r'([a-zA-Z]+)mer[^a-z]')

    def parse(self, field):
        """
        Parse uniprot field.

        Parameters
        ----------
        field : str
            field to parse.

        Returns
        -------
        int
            Stoichiometry parsed (None if field was ambiguous).

        """
        # if field is empty, stoichiometry is one
        if pandas.isnull(field):
            return None

        prefixes = self._subunit_reader.findall(field)
        # if there is only one word of the form [prefix]mer,
        if len(prefixes) == 1:
            prefix = prefixes[0].upper()
            # check for prefix in prefix rule
            return self.prefix_rule.get(prefix, None)
        else:
            return None


class CofactorParser(object):
    """Class parsing Cofactor uniprot field."""

    _name_reader = re.compile(r'Name=([^;]+); Xref=ChEBI:([^;]+);')
    _note_reader = re.compile(r'Note=(.*)')
    _stoichiometry_reader = re.compile(r'Binds ([\w]+)')

    def parse(self, field):
        """
        Parse uniprot field.

        Parameters
        ----------
        field : str
            Uniprot field containing cofactor information.

        Returns
        -------
        cofactors: list
            Cofactor objects containing information retrieved,
            where info was unambiguous.
        cofactors_to_cure : list
            Cofactor objects where some info was
            ambiguous. If some information could not be retrieved,
            its field is set to None

        """
        if pandas.isnull(field):
            return [], []
        cofactor_notes = field.split('COFACTOR:')[1:]
        cofactors = []
        needs_curation = False
        for note in cofactor_notes:
            # read name(s) and chebi identifier(s) of cofactor(s)
            # if no name was found, indicate chebi and name as missing
            full_name = self._name_reader.findall(note)
            if not full_name:
                full_name.append([None, None])
            # extract subnote if possible
            subnotes = self._note_reader.findall(note)
            subnote = subnotes[0] if len(subnotes) == 1 else note
            # infer stoichiometry:
            #  - nothing read: stoichiometry is implicitly 1
            #  - one value read: use value if can be cast to integer, else
            #    tag as missing information.
            #  - multiple values read: tag as missing information.
            stoichiometry = self._stoichiometry_reader.findall(note)
            if not stoichiometry:
                stoichiometry = 1
            elif len(stoichiometry) == 1:
                try:
                    stoichiometry = int(stoichiometry[0])
                except ValueError:
                    stoichiometry = None
            else:
                stoichiometry = None
            needs_curation = (
                needs_curation or
                (stoichiometry is None
                 or len(full_name) > 1
                 or full_name[0][0] is None)
                )
            # if there are several names, assume stoichiometry
            # is number found earlier for first element of the list
            # and 0 for the rest
            for name, chebi in full_name:
                cofactors.append(Cofactor(chebi, name, stoichiometry, subnote))
                stoichiometry = 0
        return cofactors, needs_curation
