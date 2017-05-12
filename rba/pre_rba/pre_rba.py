# global imports
import os.path

from pipeline_parameters import *
from sbml_filter import *
from uniprot_filter import *
from uniprot_importer import *
from metabolite_data import *
from macrocomponents import *
from fasta_parser import *

import rba_data
from default_processes import *
from ..rba_xml import *
from ..rba_model import *

class PreRba(object):
    def __init__(self, parameter_file):
        ## empty RBA XML structures
        self.model = RbaModel()

        ## initialize pipeline
        parameters = PipelineParameters(parameter_file).parameters
        self._organism_id = parameters['ORGANISM_ID']
        sbml_file = parameters['SBML_FILE']
        self._input_dir = parameters['INPUT_DIR']
        self.model.output_dir = parameters['OUTPUT_DIR']
        try:
            external_compartments = map(str.strip, parameters['EXTERNAL_COMPARTMENTS'].split(','))
        except KeyError:
            external_compartments = []

        ## read input data
        self._read_sbml(sbml_file, external_compartments)
        proteins_to_retrieve = [p for e in self._enzyme_comp for p in e if p != '']
        self._read_uniprot(proteins_to_retrieve)
        self.metabolite_map \
            = MetaboliteData(self._species_ids, self._cofactors,
                             self._input_dir).metabolites
        self._read_trnas()
        self._ribosome_composition = self._read_machinery('ribosome.fasta')
        self._chaperone_composition = self._read_machinery('chaperones.fasta')
        self._macro_fluxes \
            = Macrocomponents(self._species_ids,
                              self._input_dir).target_flux

        ## post processing
        self._initialize_enzymes()
        self._initialize_parameters()
        self._initialize_processes()
        # add average RNA
        rna_comp = {'A': 0.2818, 'C': 0.2181, 'G': 0.2171, 'U': 0.283}
        self.model.rnas.macromolecules.append \
            (Macromolecule(rba_data.mrna, self._cytoplasm_id, rna_comp))
        # add average DNA
        dna_comp = {'A': 0.2818, 'C': 0.2181, 'G': 0.2171, 'T': 0.283}
        self.model.dna.components.append \
            (Component('A', 'Adenosine residue', 'Nucleotide', 0))
        self.model.dna.components.append \
            (Component('C', 'Cytosine residue', 'Nucleotide', 0))
        self.model.dna.components.append \
            (Component('G', 'Guanine residue', 'Nucleotide', 0))
        self.model.dna.components.append \
            (Component('T', 'Thymine residue', 'Nucleotide', 0))
        self.model.dna.macromolecules.append(Macromolecule(rba_data.dna, self._cytoplasm_id, dna_comp))
        # medium concentrations
        for m in self._external_metabolites:        
            # /!\ we identify metabolites by their prefix !!!
            # the idea is that M_glc_p and M_glc_e will be seen
            # as the same metabolite M_glc
            self.model.medium[m.rsplit('_',1)[0]] \
                = rba_data.default_medium_concentration
            
    def _read_sbml(self, file_name, external_ids = []):
        # parse data
        sbml_data = SBMLFilter(self._input(file_name),
                               external_ids = external_ids)
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
        self.model.metabolism.compartments = ListOfCompartments()
        for c in uniprot_data.compartment_ids:
            self.model.metabolism.compartments.append(Compartment(c))
        self.model.proteins.components = uniprot_data.components
        self._cofactors = [c for c in uniprot_data.components
                           if c.type=='cofactor']
        self.model.proteins.macromolecules = uniprot_data.proteins
        self._protein_stoichiometry = uniprot_data.protein_stoichiometry
        self._cytoplasm_id = uniprot_data.cytoplasm_compartment
        self._external_id = uniprot_data.secreted_compartment
        self._average_protein_length = uniprot_data.average_protein_length
        self._sbml_to_uniprot = uniprot_data.unknown_map

    def _read_trnas(self):
        # read all real trnas (as described in fasta files)
        trna_data = FastaParser(self._input('trnas.fasta')).entries
        self.model.rnas.components.append \
            (Component('A', 'Adenosine residue', 'Nucleotide', 2.9036))
        self.model.rnas.components.append \
            (Component('C', 'Cytosine residue', 'Nucleotide', 2.7017))
        self.model.rnas.components.append \
            (Component('G', 'Guanine residue', 'Nucleotide', 3.0382))
        self.model.rnas.components.append \
            (Component('U', 'Uramine residue', 'Nucleotide', 2.7102))
        # map real trnas to user trnas
        # for example, user may agregate all trnas into a single metabolite
        # in this case, we take an average composition for a trna
        user_rnas = {}
        for rna in trna_data:
            try:
                id_ = self.metabolite_map[rna.id.upper()].sbml_id
            except KeyError:
                continue
            if not(id_): continue
            if not(user_rnas.has_key(id_)): user_rnas[id_] = []
            user_rnas[id_].append(rna.sequence)
        # write user rnas to xml structure
        for id_, seq in user_rnas.iteritems():
            comp = self._rna_composition(''.join(seq))
            average_comp = { k: float(v)/len(seq) for k, v in comp.iteritems() }
            self.model.rnas.macromolecules.append \
                (Macromolecule(id_, self._cytoplasm_id, average_comp))

    def _read_machinery(self, file_name):
        data = FastaParser(self._input(file_name)).entries
        machinery_composition = {}
        for molecule in data:
            id_ = molecule.id
            machinery_composition[id_] = molecule.stoichiometry
            if molecule.set_name == 'rna':
                comp = self._rna_composition(molecule.sequence)
                self.model.rnas.macromolecules.append \
                    (Macromolecule(id_, self._cytoplasm_id, comp))
            elif molecule.set_name == 'protein':
                comp = self._protein_composition(molecule.sequence)
                self.model.proteins.macromolecules.append \
                    (Macromolecule(id_, self._cytoplasm_id, comp))
            else:
                UserWarning(file_name + ': Unknown set ' + molecule.set_name)
        return machinery_composition
        
    def _input(self, file_name):
        return os.path.join(self._input_dir, file_name)

    def _rna_composition(self, sequence):
        comp = dict.fromkeys(['A','C','G','U'], 0)
        for n in sequence:
            if n == 'T': comp['U'] +=1
            else: comp[n] += 1
        return comp

    def _protein_composition(self, sequence):
        comp = dict.fromkeys(rba_data.aas, 0)
        for aa in sequence: comp[aa] += 1
        return comp
    
    def _initialize_enzymes(self):
        ## efficiency functions
        # for the moment simply put default value
        default_function = Function('default', 'constant', {})
        self.model.enzymes.efficiency_functions.append(default_function)

        ## enzymes
        enzyme_list = self.model.enzymes.enzymes
        def_efficiency = EnzymeEfficiency('default', {'CONSTANT': rba_data.default_catalytic_activity})
        tra_efficiency = EnzymeEfficiency('default', {'CONSTANT': rba_data.default_transporter_activity})
        MM_parameters = {'Km': rba_data.default_import_Km,
                         'kmax': rba_data.default_import_kmax}
        for reaction, comp in zip(self._reaction_ids, self._enzyme_comp):
            # base information
            enzyme = Enzyme(reaction)
            enzyme_list.append(enzyme)
            
            # machinery composition
            reactants = enzyme.machinery_composition.reactants
            for protein in comp:
                [name, sto] = self._protein_info(protein)
                if sto > 0:
                    reactants.append(SpeciesReference(name, sto))

            # base enzymatic activity
            if self._has_membrane_enzyme[reaction]:
                enzyme.enzymatic_activity.enzyme_efficiencies\
                                         .append(tra_efficiency)
            else:
                enzyme.enzymatic_activity.enzyme_efficiencies\
                                         .append(def_efficiency)

            # transport activity
            if self._imported_metabolites.has_key(reaction):
                for m in self._imported_metabolites[reaction]:
                    enzyme.enzymatic_activity.transporter_efficiency.append\
                        (Function('', 'michaelisMenten', MM_parameters, m))

        return enzyme_list

    def _enzyme_efficiency(self, function, parameters):
        new_node = self._new_node('enzymeEfficiency')
        new_node.setAttribute('function', function)
        new_node.appendChild(self._parameter_list(parameters))
        return new_node
        
    def _protein_info(self, protein):
        """
        Gather uniprot info for given protein. If no info is found, info
        for an average protein is returned.
        """
        # check in protein name has been overriden
        try:
            protein = self._sbml_to_uniprot[protein]
        except KeyError:
            pass
        # check if reaction is spontaneous
        if protein == '': return [protein, 0]
        # check if protein is known
        try:
            return [protein, self._protein_stoichiometry[protein]]
        except KeyError:
            # if protein is unknown we replace it with an average protein
            return [rba_data.average_protein_id(self._cytoplasm_id), 1]

    def _initialize_parameters(self):
        ## maximal densities
        other_compartments = self._compartment_ids[:]
        other_compartments.remove(self._cytoplasm_id)
        other_compartments.remove(self._external_id)
        constraints = self.model.parameters.maximal_densities
        # always start with cytoplasm
        c_density = MaximalDensity(self._cytoplasm_id)
        c_density.value = rba_data.cytoplasm_density
        constraints.append(c_density)
        # now treat the rest
        for user_id in other_compartments:
            new_density = MaximalDensity(user_id)
            new_density.function_references.append('amino_acid_concentration')
            new_density.function_references.append(rba_data.protein_fraction_id(user_id))
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
        for fn in rba_data.fns: fns.append(fn)

    def _initialize_processes(self):
        ## create processes
        self.model.processes = DefaultProcesses \
                         (self._ribosome_composition,
                          self._chaperone_composition,
                          self._compartment_ids, self._cofactors,
                          self.metabolite_map, self._macro_fluxes) \
                          .rba_processes
        ## add atpm reaction
        reaction = Reaction(rba_data.atpm_reaction, False)
        for m in ['ATP', 'H2O']:
            id_ = self.metabolite_map[m].sbml_id
            if id_: reaction.reactants.append(SpeciesReference(id_, 1))
        for m in ['ADP', 'H', 'Pi']:
            id_ = self.metabolite_map[m].sbml_id
            if id_: reaction.products.append(SpeciesReference(id_, 1))
        self.model.metabolism.reactions.append(reaction)
        enzyme = Enzyme(rba_data.atpm_reaction)
        enzyme.enzymatic_activity.enzyme_efficiencies.append(EnzymeEfficiency('default', {'CONSTANT': rba_data.default_catalytic_activity}))
        self.model.enzymes.enzymes.append(enzyme)
                
