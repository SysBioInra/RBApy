"""
Module defining PreRba class.
"""

# python 2/3 compatibility
from __future__ import absolute_import, division, print_function

# global imports
import os.path

# data processing imports
from rba.pre_rba.pipeline_parameters import PipelineParameters
from rba.pre_rba.sbml_filter import SBMLFilter
from rba.pre_rba.uniprot_filter import UniprotFilter
from rba.pre_rba.uniprot_importer import UniprotImporter
from rba.pre_rba.metabolite_data import MetaboliteData
from rba.pre_rba.macrocomponents import Macrocomponents
from rba.pre_rba.fasta_parser import FastaParser

# default data imports
from rba.pre_rba import rba_data
from rba.pre_rba.default_processes import DefaultProcesses

# model imports
from rba import rba_xml
from rba.rba_model import RbaModel

class PreRba(object):
    """
    Class generating RBA model from SBML file and helper files.

    Attributes:
        model: RbaModel generated.
        metabolite_map: dictionary mapping standard metabolite names
            (as defined in rba_dat) to sbml identifiers.
    """

    def __init__(self, parameter_file):
        """
        Constructor from parameter_file.

        Args:
            parameter_file: standard pipeline input file (containing e.g.
                sbml location, organism id).
        """
        ## empty RBA XML structures
        self.model = RbaModel()

        ## initialize pipeline
        parameters = PipelineParameters(parameter_file).parameters
        self._organism_id = parameters['ORGANISM_ID']
        sbml_file = parameters['SBML_FILE']
        self._input_dir = parameters['INPUT_DIR']
        self.model.output_dir = parameters['OUTPUT_DIR']
        external_line = parameters.get('EXTERNAL_COMPARTMENTS', None)
        if external_line is not None:
            external_compartments = external_line.split(',')
            external_compartments = [e.strip() for e in external_compartments]
        else:
            external_compartments = []

        ## read input data
        self._read_sbml(sbml_file, external_compartments)
        proteins_to_retrieve = []
        for enzyme in self._enzyme_comp:
            proteins_to_retrieve += [p for p in enzyme if p != '']
        self._read_uniprot(proteins_to_retrieve)
        self.metabolite_map \
            = MetaboliteData(self._species_ids, self._cofactors,
                             self._input_dir).metabolites
        self._read_trnas()
        self._ribosome_composition = self._read_machinery('ribosome.fasta')
        self._chaperone_composition = self._read_machinery('chaperones.fasta')
        self._macro_fluxes \
            = Macrocomponents(self._species_ids, self._input_dir).target_flux

        ## post processing
        self._initialize_enzymes()
        self._initialize_parameters()
        self._initialize_processes()
        # add average RNA
        rna_comp = {'A': 0.2818, 'C': 0.2181, 'G': 0.2171, 'U': 0.283}
        rna = rba_xml.Macromolecule(rba_data.mrna, self._cytoplasm_id, rna_comp)
        self.model.rnas.macromolecules.append(rna)
        # add average DNA
        dna_comp = {'A': 0.2818, 'C': 0.2181, 'G': 0.2171, 'T': 0.283}
        comp_a = rba_xml.Component('A', 'Adenosine residue', 'Nucleotide', 0)
        comp_c = rba_xml.Component('C', 'Cytosine residue', 'Nucleotide', 0)
        comp_g = rba_xml.Component('G', 'Guanine residue', 'Nucleotide', 0)
        comp_t = rba_xml.Component('T', 'Thymine residue', 'Nucleotide', 0)
        for comp in [comp_a, comp_c, comp_g, comp_t]:
            self.model.dna.components.append(comp)
        dna = rba_xml.Macromolecule(rba_data.dna, self._cytoplasm_id, dna_comp)
        self.model.dna.macromolecules.append(dna)
        # medium concentrations
        for metab in self._external_metabolites:
            # !!! we identify metabolites by their prefix !!!
            # the idea is that M_glc_p and M_glc_e will be seen
            # as the same metabolite M_glc
            self.model.medium[metab.rsplit('_', 1)[0]] \
                = rba_data.default_medium_concentration

    def _read_sbml(self, file_name, external_ids=None):
        """
        Read sbml file and store relevant information.

        Args:
            file_name: name of sbml file.
            external_ids: sbml ids of external compartments.
        """
        if external_ids is None:
            external_ids = []
        # parse data
        sbml_data = SBMLFilter(self._input(file_name),
                               external_ids=external_ids)
        # store data
        self.model.metabolism.species = sbml_data.species
        self.model.metabolism.reactions = sbml_data.reactions
        self._enzyme_comp = sbml_data.enzymes
        self._external_metabolites = sbml_data.external_metabolites
        self._imported_metabolites = sbml_data.imported_metabolites
        self._has_membrane_enzyme = sbml_data.has_membrane_enzyme
        self._species_ids = [s.id for s in sbml_data.species]
        self._reaction_ids = [r.id for r in sbml_data.reactions]

    def _read_uniprot(self, protein_ids):
        """
        Read uniprot file and store relevant information.

        Args:
            protein_ids: identifiers of proteins to retrieve in uniprot.
        """
        input_file = self._input('uniprot.csv')
        # retrieve uniprot data from the internet if not provided by the user
        if os.path.isfile(input_file):
            print('Found uniprot file...')
        else:
            print('Could not find uniprot file. Downloading most recent'
                  ' version...')
            raw_data = UniprotImporter(self._organism_id).data
            if len(raw_data) == 0:
                raise UserWarning('Invalid organism, could not retrieve Uniprot'
                                  ' data. Interrupting execution.')
            with open(input_file, 'w') as f:
                f.write(raw_data)
        # parse data
        uniprot_data = UniprotFilter(protein_ids, self._input_dir)

        # store data
        self._compartment_ids = uniprot_data.compartment_ids
        self.model.metabolism.compartments = rba_xml.ListOfCompartments()
        for c in uniprot_data.compartment_ids:
            self.model.metabolism.compartments.append(rba_xml.Compartment(c))
        self.model.proteins.components = uniprot_data.components
        self._cofactors = [c for c in uniprot_data.components
                           if c.type == 'cofactor']
        self.model.proteins.macromolecules = uniprot_data.proteins
        self._protein_stoichiometry = uniprot_data.protein_stoichiometry
        self._cytoplasm_id = uniprot_data.cytoplasm_compartment
        self._external_id = uniprot_data.secreted_compartment
        self._average_protein_length = uniprot_data.average_protein_length
        self._sbml_to_uniprot = uniprot_data.unknown_map

    def _read_trnas(self):
        """
        Read trna in fasta file and store them in the model.
        """
        # read all real trnas (as described in fasta files)
        trna_data = FastaParser(self._input('trnas.fasta')).entries
        A = rba_xml.Component('A', 'Adenosine residue', 'Nucleotide', 2.9036)
        C = rba_xml.Component('C', 'Cytosine residue', 'Nucleotide', 2.7017)
        G = rba_xml.Component('G', 'Guanine residue', 'Nucleotide', 3.0382)
        U = rba_xml.Component('U', 'Uramine residue', 'Nucleotide', 2.7102)
        for comp in [A, C, G, U]:
            self.model.rnas.components.append(comp)
        # map real trnas to user trnas
        # for example, user may agregate all trnas into a single metabolite
        # in this case, we take an average composition for a trna
        user_rnas = {}
        for rna in trna_data:
            metab = self.metabolite_map.get(rna.id.upper(), None)
            if metab is not None:
                id_ = metab.sbml_id
                if id_:
                    rna_list = user_rnas.setdefault(id_, [])
                    rna_list.append(rna.sequence)
        # write user rnas to xml structure
        for id_, seq in user_rnas.items():
            comp = self._rna_composition(''.join(seq))
            average_comp = {k: float(v)/len(seq) for k, v in comp.items()}
            rna = rba_xml.Macromolecule(id_, self._cytoplasm_id, average_comp)
            self.model.rnas.macromolecules.append(rna)

    def _read_machinery(self, file_name):
        """
        Read machinery in fasta file and return its composition.

        Returns:
            A dict mapping components of machinery to their stoichiometry.
        """
        data = FastaParser(self._input(file_name)).entries
        machinery_composition = {}
        for molecule in data:
            id_ = molecule.id
            machinery_composition[id_] = molecule.stoichiometry
            if molecule.set_name == 'rna':
                comp = self._rna_composition(molecule.sequence)
                rna = rba_xml.Macromolecule(id_, self._cytoplasm_id, comp)
                self.model.rnas.macromolecules.append(rna)
            elif molecule.set_name == 'protein':
                comp = self._protein_composition(molecule.sequence)
                prot = rba_xml.Macromolecule(id_, self._cytoplasm_id, comp)
                self.model.proteins.macromolecules.append(prot)
            else:
                UserWarning(file_name + ': Unknown set ' + molecule.set_name)
        return machinery_composition

    def _input(self, file_name):
        """
        Return full path to file located in input directory.
        """
        return os.path.join(self._input_dir, file_name)

    @staticmethod
    def _rna_composition(sequence):
        """
        Compute composition of rna sequence.
        """
        comp = dict.fromkeys(['A', 'C', 'G', 'U'], 0)
        for n in sequence:
            if n == 'T':
                comp['U'] += 1
            else:
                comp[n] += 1
        return comp

    @staticmethod
    def _protein_composition(sequence):
        """
        Compute composition of protein sequence.
        """
        comp = dict.fromkeys(rba_data.aas, 0)
        for aa in sequence:
            comp[aa] += 1
        return comp

    def _initialize_enzymes(self):
        """
        Create enzyme sturcture of RBA model.
        """
        ## efficiency functions
        # for the moment simply put default value
        default_function = rba_xml.Function('default', 'constant', {})
        self.model.enzymes.efficiency_functions.append(default_function)

        ## enzymes
        enzyme_list = self.model.enzymes.enzymes
        def_param = {'CONSTANT': rba_data.default_catalytic_activity}
        def_efficiency = rba_xml.EnzymeEfficiency('default', def_param)
        tra_param = {'CONSTANT': rba_data.default_transporter_activity}
        tra_efficiency = rba_xml.EnzymeEfficiency('default', tra_param)
        mm_parameters = {'Km': rba_data.default_import_Km,
                         'kmax': rba_data.default_import_kmax}
        for reaction, comp in zip(self._reaction_ids, self._enzyme_comp):
            # base information
            enzyme = rba_xml.Enzyme(reaction)
            enzyme_list.append(enzyme)
            # machinery composition
            reactants = enzyme.machinery_composition.reactants
            for protein in comp:
                [name, sto] = self._protein_info(protein)
                if sto > 0:
                    reactants.append(rba_xml.SpeciesReference(name, sto))
            # base enzymatic activity
            efficiencies = enzyme.enzymatic_activity.enzyme_efficiencies
            if self._has_membrane_enzyme[reaction]:
                efficiencies.append(tra_efficiency)
            else:
                efficiencies.append(def_efficiency)
            # transport activity
            transport = enzyme.enzymatic_activity.transporter_efficiency
            if reaction in self._imported_metabolites:
                for m in self._imported_metabolites[reaction]:
                    transport.append(rba_xml.Function('', 'michaelisMenten',
                                                      mm_parameters, m))
        return enzyme_list

    def _protein_info(self, protein):
        """
        Return uniprot name and stoichiometry for protein defined in SBML.

        Args:
            protein: sbml id of protein as indicated in enzyme notes.

        Returns:
            List where first element is uniprot name of protein and second
            element is stoichiometry. If no uniprot info was found, info
            for an average protein is returned.
        """
        # retrieve protein name if it has been overriden
        protein = self._sbml_to_uniprot.get(protein, protein)
        # check if reaction is spontaneous
        if not protein:
            return [protein, 0]
        # check if protein is known
        # else replace it with an average protein in the cytosol
        sto = self._protein_stoichiometry.get(protein, None)
        if sto is not None:
            return [protein, sto]
        else:
            return [rba_data.average_protein_id(self._cytoplasm_id), 1]

    def _initialize_parameters(self):
        """
        Create parameter structure of RBA model.
        """
        ## target densities
        other_compartments = self._compartment_ids[:]
        other_compartments.remove(self._cytoplasm_id)
        other_compartments.remove(self._external_id)
        constraints = self.model.parameters.target_densities
        aggregates = self.model.parameters.aggregates
        # always start with cytoplasm
        c_density = rba_xml.TargetDensity(self._cytoplasm_id)
        c_density.upper_bound = rba_data.cytoplasm_density
        constraints.append(c_density)
        # now treat the rest
        for user_id in other_compartments:
            # create aggregated function for density
            id_ = user_id + '_density'
            new_aggregate = rba_xml.Aggregate(id_, 'multiplication')
            for ref in ['amino_acid_concentration',
                        rba_data.protein_fraction_id(user_id)]:
                new_aggregate.function_references \
                             .append(rba_xml.FunctionReference(ref))
            aggregates.append(new_aggregate)
            # create constraint
            new_density = rba_xml.TargetDensity(user_id)
            new_density.upper_bound = id_
            constraints.append(new_density)

        ## functions
        fns = self.model.parameters.functions
        # protein fractions
        fns.append(rba_data.protein_fraction(self._cytoplasm_id,
                                             rba_data.cytoplasm_fraction))
        fns.append(rba_data.protein_fraction(self._external_id,
                                             rba_data.secreted_fraction))
        other_fraction = rba_data.other_fraction / len(other_compartments)
        for c in other_compartments:
            fns.append(rba_data.protein_fraction(c, other_fraction))
        # non enzymatic fractions
        fns.append(rba_data.non_enzymatic_fraction_cytoplasm(self._cytoplasm_id))
        fns.append(rba_data.non_enzymatic_fraction_secreted(self._external_id))
        for c in other_compartments:
            fns.append(rba_data.non_enzymatic_fraction_other(c))
        # protein length
        fns.append(rba_data.inverse_average_protein_length(self._average_protein_length))
        # add all other functions
        for fn in rba_data.fns:
            fns.append(fn)

    def _initialize_processes(self):
        """
        Create process structure of RBA model.
        """
        ## create processes
        processes = DefaultProcesses(self._ribosome_composition,
                                     self._chaperone_composition,
                                     self._compartment_ids, self._cofactors,
                                     self.metabolite_map, self._macro_fluxes)
        # add processes to model
        self.model.processes = processes.rba_processes
        # add aggregated functions to model
        for agg in processes.aggregates:
            self.model.parameters.aggregates.append(agg)
        ## add atpm reaction
        reaction = rba_xml.Reaction(rba_data.atpm_reaction, False)
        for m in ['ATP', 'H2O']:
            id_ = self.metabolite_map[m].sbml_id
            if id_:
                reaction.reactants.append(rba_xml.SpeciesReference(id_, 1))
        for m in ['ADP', 'H', 'Pi']:
            id_ = self.metabolite_map[m].sbml_id
            if id_:
                reaction.products.append(rba_xml.SpeciesReference(id_, 1))
        self.model.metabolism.reactions.append(reaction)
        enzyme = rba_xml.Enzyme(rba_data.atpm_reaction)
        activity = {'CONSTANT': rba_data.default_catalytic_activity}
        efficiency = rba_xml.EnzymeEfficiency('default', activity)
        enzyme.enzymatic_activity.enzyme_efficiencies.append(efficiency)
        self.model.enzymes.enzymes.append(enzyme)
