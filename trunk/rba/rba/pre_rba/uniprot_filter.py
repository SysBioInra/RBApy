
# make division python3 style (float division), i.e. 3 / 2 = 1.5
from __future__ import division

import pandas, os.path
import re

from ..rba_xml import *
from .cofactor import *
from .subunit import *
from .location import *
from .unknown_proteins import *

class UniprotFilter:
    """
    Class used to filter RBA-relevant data from a uniprot file.

    :param input_dir: Directory containing uniprot and helper files.
    :param gene_ids: Genes which information needs to be retrieved. If nothing
    is specified, everything is kept.
    :type input_dir: string
    :type gene_ids: list of strings
    """
    def __init__(self, gene_ids, input_dir = '.'):
        """
        Constructor.

        :param gene_ids: list of proteins that need to be retrieved.
        :param input_dir: path to folder containing uniprot files.
        :type gene_ids: list of strings
        :type input_dir: string
        """
        self._gene_ids = list(set(gene_ids))
        self._upper_gene_ids = list(map(str.upper, self._gene_ids))
        self._gene_id_found = [False] * len(self._gene_ids)
        
        # open uniprot data
        input_file = os.path.join(input_dir, 'uniprot.csv')
        data = pandas.read_csv(input_file, sep = '\t')
        # open helper file for unknown proteins
        unknown_proteins = UnknownProteins(input_dir)
        self.unknown_map = unknown_proteins.data
        # replace ids marked as unknown by user-define counterpart
        for i, gene in enumerate(gene_ids):
            try:
                gene_ids[i] = self.unknown_map[gene]
            except KeyError: pass
            
        # transform gene data into list
        self._gene_reader = re.compile(r'([^\s]+)')
        data['SBML name'] = data['Gene names'].apply(self._genes_as_list)
        
        # keep only ids that appear in the SBML file or helper file
        data['SBML name'] = data['SBML name'].apply(self._filter_genes)

        # update helper file for unknown proteins
        not_found = []
        for gene_id, found in zip(self._gene_ids, self._gene_id_found):
            if not(found):
                not_found.append(gene_id)
        unknown_proteins.add(not_found)

        # remove entries for which there is no known id
        drop_rows = []
        gene_data = data['SBML name']
        for i, gene_names in enumerate(gene_data):
            if len(gene_names) == 0:
                drop_rows.append(i)
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
        self.components = ListOfComponents()
        self.proteins = ListOfMacromolecules()
        self.protein_stoichiometry = {}
        # add amino acids to component set
        amino_acids = set()
        for comp in aa_data.values():
            for aa in comp:
                amino_acids.add(aa)
        for aa in amino_acids:
            self.components.append(Component(aa, '', 'amino_acid', 1))
        # store protein data
        entry_data = reduced_data['Entry']
        gene_data = reduced_data['SBML name']
        known_cofactors = []
        for entry, gene_names in zip(entry_data, gene_data):
            # add amino acids to composition
            composition = aa_data[entry].copy()
            # add cofactors to composition
            cofactors = []
            if entry in cofactor_data:
                cofactors = cofactor_data[entry]
            for c in cofactors:
                composition[c.chebi] = c.stoichiometry
                if c.chebi not in known_cofactors:
                    self.components.append(Component(c.chebi, c.name,
                                                     'cofactor', 0))
                    known_cofactors.append(c.chebi)
            # check if protein is referenced in SBML with different ids (!!!)
            if len(gene_names) > 1:
                print('In your SBML, ids ' + ', '.join(gene_names)
                      + ' reference the same protein. Protein will be '
                      + 'duplicated for now. Check your SBML for future runs.')
            for gene_name in gene_names:
                self.proteins.append \
                    (Macromolecule(gene_name, location_data[entry], composition))
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
            self.proteins.append \
                (Macromolecule(name, location, average_protein))
            self.protein_stoichiometry[name] = 1

    def _genes_as_list(self, gene_names):
        """
        Transforms gene data in a format that is easier to work on.
        """
        # no gene name associated with protein
        if pandas.isnull(gene_names): return []

        # transform names into a list
        return list(set(map(str.upper, self._gene_reader.findall(gene_names))))

    def _filter_genes(self, gene_ids):
        # filter known ids
        # we expect to keep at most one !!!
        curated_ids = []
        for i in gene_ids:
            try:
                index = self._upper_gene_ids.index(i)
                curated_ids.append(self._gene_ids[index])
                self._gene_id_found[index] = True
            except ValueError:
                pass
        return curated_ids

    def _aa_data(self, uniprot_data):
        entry_data = uniprot_data['Entry']
        sequence_data = uniprot_data['Sequence']
        aa_data = {}
        # read sequence composition for every entry
        for entry, sequence in zip(entry_data, sequence_data):
            new_data = {}
            for c in sequence:
                if c not in new_data:
                    new_data[c] = 0
                new_data[c] += 1
            aa_data[entry] = new_data
        return aa_data
