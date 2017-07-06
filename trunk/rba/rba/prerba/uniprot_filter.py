"""
Module defining UniprotFilter class.
"""

# python 2/3 compatibility
from __future__ import division, print_function

# global imports
import os.path
import re
import pandas

# local imports
import rba.xml
from rba.prerba.cofactor import CofactorData
from rba.prerba.subunit import SubunitData
from rba.prerba.location import LocationData
from rba.prerba.unknown_proteins import UnknownProteins

class UniprotFilter(object):
    """
    Class used to filter RBA-relevant data from a uniprot file.

    Attributes:
        unknown_map: dict mapping sbml genes not retrieved in uniprot with
            a corresponding user-defined gene theoretically found in uniprot.
        default_location: name of default compartment.
        cytoplasm_compartment: name of cytoplasm.
        secreted_compartment: name of external compartment.
        compartment_ids: list of names of all compartments.
        components: rba.xml.ListOfComponents containing all elements that are
            used to build proteins (amino acids, vitamins, cofactors).
        proteins: rba.xml.ListOfMacromolecules containing all proteins that
            have been filtered.
        protein_stoichiometry: dict mapping protein identifiers with their
            stoichiometry within their enzymatic complex.
        average_protein_length: average protein length.
    """

    def __init__(self, gene_ids, input_dir='.'):
        """
        Constructor from input directory and gene identifiers.

        Args:
            input_dir: Directory containing uniprot and helper files.
            gene_ids: Genes to retrieve. If nothing is specified, everything
                is kept.
        """
        # open uniprot data
        input_file = os.path.join(input_dir, 'uniprot.csv')
        data = pandas.read_csv(input_file, sep='\t')
        # add column containing identifier of protein as found in gene_ids
        unknown_proteins = UnknownProteins(input_dir)
        self.unknown_map = unknown_proteins.data
        data['Lookup names'], not_found \
            = self._lookup_names(data['Gene names'], gene_ids, self.unknown_map)
        unknown_proteins.add(not_found)

        # remove entries for which there is no known id
        drop_rows = []
        gene_data = data['Lookup names']
        for index, gene_names in enumerate(gene_data):
            if len(gene_names) == 0:
                drop_rows.append(index)
        reduced_data = data.drop(drop_rows)

        # get cofactor, subunit and location data
        cofactor_data = CofactorData(reduced_data, input_dir).cofactors
        subunit_data = SubunitData(reduced_data, input_dir).subunits
        location_parser = LocationData(reduced_data, input_dir)
        location_data = location_parser.locations
        self.default_location = location_parser.default_location
        self.cytoplasm_compartment = location_parser.cytoplasm_id()
        self.secreted_compartment = location_parser.secreted_id()
        self.compartment_ids = location_parser.compartment_ids()

        # gather protein compositions in terms of amino acids
        aa_data = self._aa_data(data)

        # store data
        self.components = rba.xml.ListOfComponents()
        self.proteins = rba.xml.ListOfMacromolecules()
        self.protein_stoichiometry = {}
        # add amino acids to component set
        amino_acids = set()
        for comp in aa_data.values():
            for aa in comp:
                amino_acids.add(aa)
        for aa in amino_acids:
            self.components.append(rba.xml.Component(aa, '', 'amino_acid', 1))
        # store protein data
        entry_data = reduced_data['Entry']
        gene_data = reduced_data['Lookup names']
        known_cofactors = []
        for entry, gene_names in zip(entry_data, gene_data):
            # add amino acids to composition
            composition = aa_data[entry].copy()
            # add cofactors to composition
            cofactors = []
            if entry in cofactor_data:
                cofactors = cofactor_data[entry]
            for cofactor in cofactors:
                composition[cofactor.chebi] = cofactor.stoichiometry
                if cofactor.chebi not in known_cofactors:
                    new_cofactor = rba.xml.Component(cofactor.chebi,
                                                     cofactor.name,
                                                     'cofactor', 0)
                    self.components.append(new_cofactor)
                    known_cofactors.append(cofactor.chebi)
            # check if protein is referenced in SBML with different ids (!!!)
            if len(gene_names) > 1:
                print('In your SBML, ids ' + ', '.join(gene_names)
                      + ' reference the same protein. Protein will be '
                      + 'duplicated for now. Check your SBML for future runs.')
            # add protein(s) to list
            for gene_name in gene_names:
                protein = rba.xml.Macromolecule(gene_name, location_data[entry],
                                                composition)
                self.proteins.append(protein)
                self.protein_stoichiometry[gene_name] = subunit_data[entry]

        # compute average protein
        average_protein = dict.fromkeys(amino_acids, 0)
        for composition in aa_data.values():
            for aa in composition:
                if aa in average_protein:
                    average_protein[aa] += composition[aa]
        for aa in average_protein:
            average_protein[aa] /= len(aa_data)
        self.average_protein_length = sum(average_protein.values())
        # add one average protein for every compartment
        for location in self.compartment_ids:
            name = 'average_protein_' + location
            protein = rba.xml.Macromolecule(name, location, average_protein)
            self.proteins.append(protein)
            self.protein_stoichiometry[name] = 1

    @staticmethod
    def _lookup_names(uniprot_names, gene_ids, unknown_map):
        """
        Create column containing lookup name of proteins.

        Args:
            uniprot_names: column containing raw uniprot 'Gene names' fields.
            gene_ids: name of genes to retrieve.
            unknown_map: dict where keys are genes that could not be 
                retrieved in 'Gene names'. Values is an alternative identifier
                which might be found in 'Gene names'.

        Returns:
            First element: Column containing list of genes found on every
                row of input column, using names in gene_ids.
                Example:
                    uniprot_names = [['gene_a gene_a_2'],
                                     ['gene_b]]
                    gene_ids = ['GENE_A', 'GENE_C']
                    output -> [['GENE_A'],
                               []]
            Second element: List of ids that could not be retrieved.
        """
        gene_ids = list(set(gene_ids))
        ## replace ids by user-defined counterpart if applicable
        for index, name in enumerate(gene_ids):
            user_id = unknown_map.get(name, None)
            if user_id is not None:
                gene_ids[index] = user_id
        ## create a mapping with standardized upper case names
        original_gene_ids = {name.upper(): name for name in gene_ids}
        ## match uniprot field with genes to retrieve
        gene_reader = re.compile(r'([^\s]+)')
        lookup_names = uniprot_names.copy()
        gene_found = {name: False for name in gene_ids}
        for row_index, row in enumerate(lookup_names):
            # transform raw uniprot field into standardized list
            if pandas.isnull(row):
                gene_list = []
            else:
                gene_list = list(set(gene_reader.findall(row)))
            # keep only ids that appear in the list of ids to retrieve
            curated_ids = []
            for name in gene_list:
                original_name = original_gene_ids.get(name.upper(), None)
                if original_name is not None:
                    curated_ids.append(original_name)
                    gene_found[original_name] = True
            lookup_names[row_index] = curated_ids
        ## establish list of names that could not be found
        not_found = [name for name, found in gene_found.items() if not found]
        return lookup_names, not_found

    @staticmethod
    def _aa_data(uniprot_data):
        """
        Retrive amino acid composition of proteins.

        Args:
            uniprot_data: pandas table containing protein information.

        Returns:
            dict mapping uniprot entry to its composition. Composition is a
            dict mapping amino acids with their occurence in the sequence.
        """
        entry_data = uniprot_data['Entry']
        sequence_data = uniprot_data['Sequence']
        aa_data = {}
        # read sequence composition for every entry
        for entry, sequence in zip(entry_data, sequence_data):
            new_data = {}
            for char in sequence:
                new_data.setdefault(char, 0)
                new_data[char] += 1
            aa_data[entry] = new_data
        return aa_data
