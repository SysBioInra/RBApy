"""Model used to build RBA XML model."""

# python 2/3 compatibility
from __future__ import division, print_function, absolute_import

# global imports
import itertools
import copy

# local imports
from rba.prerba.user_data import ntp_composition
from rba.prerba.default_processes import DefaultProcesses
import rba.xml


class ModelBuilder(object):

    def __init__(self, default_data, user_data):
        self.data = user_data
        self.default = default_data

    def build_metabolism(self):
        """
        Build metabolism part of RBA model.

        Returns
        -------
        rba.xml.RbaMetabolism
            RBA metabolism model in XML format.
        """
        metabolism = rba.xml.RbaMetabolism()

        metabolism.species = copy.deepcopy(self.data.sbml_species())
        metabolism.reactions = copy.deepcopy(self.data.sbml_reactions())
        # add atpm reaction
        reaction = rba.xml.Reaction(self.default.atpm_reaction, False)
        for m in ['ATP', 'H2O']:
            id_ = self.data.metabolite_map[m].sbml_id
            if id_:
                reaction.reactants.append(rba.xml.SpeciesReference(id_, 1))
        for m in ['ADP', 'H', 'Pi']:
            id_ = self.data.metabolite_map[m].sbml_id
            if id_:
                reaction.products.append(rba.xml.SpeciesReference(id_, 1))
        metabolism.reactions.append(reaction)
        # add compartments
        for c in self.data.compartments():
            metabolism.compartments.append(rba.xml.Compartment(c))
        return metabolism

    def build_density(self):
        """
        Build density part of RBA model.

        Returns
        -------
        rba.xml.RbaDensity
            RBA density model in XML format.

        """
        density = rba.xml.RbaDensity()
        constraints = density.target_densities
        external = self.data.compartment('Secreted')
        for c_id in self.data.compartments():
            if c_id != external:
                new_density = rba.xml.TargetDensity(c_id)
                new_density.upper_bound = c_id + '_density'
                constraints.append(new_density)
        return density

    def build_parameters(self):
        """
        Build parameter part of RBA model.

        Returns
        -------
        rba.xml.RbaParameters
            RBA parameter model in XML format.

        """
        parameters = rba.xml.RbaParameters()
        def_params = self.default.parameters
        # density related functions
        cytoplasm = self.data.compartment('Cytoplasm')
        external = self.data.compartment('Secreted')
        other_cpt = self.data.compartments()
        other_cpt.remove(cytoplasm)
        other_cpt.remove(external)
        fns, aggs = def_params.density_functions(cytoplasm, external,
                                                 other_cpt)
        for fn in fns:
            parameters.functions.append(fn)
        for agg in aggs:
            parameters.aggregates.append(agg)
        # protein length
        parameters.functions.append(
            def_params.inverse_average_protein_length(
                sum(self.data.average_protein.values())
                )
            )
        # process functions
        fns, aggs = def_params.process_functions()
        for fn in fns:
            parameters.functions.append(fn)
        for agg in aggs:
            parameters.aggregates.append(agg)
        # target functions for metabolites and macrocomponents
        for id_, concentration in self.data.macrocomponents:
            parameters.functions.append(
                def_params.metabolite_concentration_function(
                    id_, concentration
                    )
                )
        for metabolite in self.data.metabolite_map.values():
            if metabolite.sbml_id and metabolite.concentration:
                parameters.functions.append(
                    def_params.metabolite_concentration_function(
                        metabolite.sbml_id, metabolite.concentration
                        )
                    )
        # enzyme efficiencies
        # base enzymatic activity
        parameters.functions.append(
            self.default.activity.efficiency_function()
            )
        parameters.functions.append(
            self.default.activity.transport_function()
            )
        # transport functions and aggregates
        reaction_ids = [r.id for r in self.data.sbml_reactions()]
        for reaction in reaction_ids:
            if self.data.has_membrane_enzyme(reaction):
                fns, agg = self.default.activity.transport_aggregate(
                    reaction, self.data.imported_metabolites(reaction)
                    )
                for fn in fns:
                    parameters.functions.append(fn)
                parameters.aggregates.append(agg)
        return parameters

    def build_proteins(self):
        """
        Build protein part of RBA model.

        Returns
        -------
        rba.xml.RbaMacromolecules
            RBA protein model in XML format.

        """
        proteins = rba.xml.RbaMacromolecules()
        # components
        for aa in self.default.metabolites.aas:
            proteins.components.append(
                rba.xml.Component(aa, '', 'amino_acid', 1)
                )
        for c in self.data.cofactors():
            proteins.components.append(
                rba.xml.Component(c.chebi, c.name, 'cofactor', 0)
                )
        # enzymatic proteins
        for gene_name, protein in self.data.enzymatic_proteins.items():
            comp = self.data.aa_composition(protein.sequence)
            for cofactor in protein.cofactors:
                comp[cofactor.chebi] = cofactor.stoichiometry
            proteins.macromolecules.append(
                rba.xml.Macromolecule(gene_name, protein.location, comp)
                )
        # average proteins
        for compartment in self.data.compartments():
            id_ = self.data.average_protein_id(compartment)
            proteins.macromolecules.append(
                rba.xml.Macromolecule(id_, compartment,
                                      self.data.average_protein)
                )
        # machinery proteins
        cytoplasm = self.data.compartment('Cytoplasm')
        for molecule in itertools.chain(self.data.ribosome,
                                        self.data.chaperone):
            if molecule.set_name == 'protein':
                proteins.macromolecules.append(
                    rba.xml.Macromolecule(
                        molecule.id, cytoplasm,
                        self.data.aa_composition(molecule.sequence)
                    )
                )
        return proteins

    def build_rnas(self):
        """
        Build RNA part of RBA model.

        Returns
        -------
        rba.xml.RbaMacromolecules
            RBA RNA model in XML format.

        """
        rnas = rba.xml.RbaMacromolecules()
        # components
        rnas.components.append(
            rba.xml.Component('A', 'Adenosine residue', 'Nucleotide', 2.9036)
            )
        rnas.components.append(
            rba.xml.Component('C', 'Cytosine residue', 'Nucleotide', 2.7017)
            )
        rnas.components.append(
            rba.xml.Component('G', 'Guanine residue', 'Nucleotide', 3.0382)
            )
        rnas.components.append(
            rba.xml.Component('U', 'Uramine residue', 'Nucleotide', 2.7102)
            )
        # user rnas
        cytoplasm = self.data.compartment('Cytoplasm')
        for rna_id, composition in self.data.rna_data.items():
            rnas.macromolecules.append(rba.xml.Macromolecule(
                rna_id, cytoplasm, composition
                ))
        # average RNA
        rnas.macromolecules.append(rba.xml.Macromolecule(
            self.default.metabolites.mrna, cytoplasm,
            {'A': 0.2818, 'C': 0.2181, 'G': 0.2171, 'U': 0.283}
            ))
        # machinery rnas
        for molecule in itertools.chain(self.data.ribosome,
                                        self.data.chaperone):
            if molecule.set_name == 'rna':
                rnas.macromolecules.append(rba.xml.Macromolecule(
                    molecule.id, cytoplasm,
                    ntp_composition(molecule.sequence)
                    ))
        return rnas

    def build_dna(self):
        """
        Build DNA part of RBA model.

        Returns
        -------
        rba.xml.RbaMacromolecules
            RBA DNA model in XML format.

        """
        dna = rba.xml.RbaMacromolecules()
        # components
        comp_a = rba.xml.Component('A', 'Adenosine residue', 'Nucleotide', 0)
        comp_c = rba.xml.Component('C', 'Cytosine residue', 'Nucleotide', 0)
        comp_g = rba.xml.Component('G', 'Guanine residue', 'Nucleotide', 0)
        comp_t = rba.xml.Component('T', 'Thymine residue', 'Nucleotide', 0)
        for comp in [comp_a, comp_c, comp_g, comp_t]:
            dna.components.append(comp)
        # average DNA
        dna.macromolecules.append(
            rba.xml.Macromolecule(
                self.default.metabolites.dna,
                self.data.compartment('Cytoplasm'),
                {'A': 0.2818, 'C': 0.2181, 'G': 0.2171, 'T': 0.283}
                )
            )
        return dna

    def build_enzymes(self):
        """
        Build enzyme part of RBA model.

        Returns
        -------
        rba.xml.RbaEnzymes
            RBA enzyme model in XML format.

        """
        enzymes = rba.xml.RbaEnzymes()
        # user enzymes
        def_act = self.default.activity
        reaction_ids = [r.id for r in self.data.sbml_reactions()]
        for reaction, comp in zip(reaction_ids,
                                  self.data.enzyme_composition()):
            if self.data.has_membrane_enzyme(reaction):
                forward = def_act.transport_aggregate_id(reaction)
                backward = def_act.transport_id
            else:
                forward = backward = def_act.efficiency_id
            # base information
            enzyme = rba.xml.Enzyme(reaction + '_enzyme', reaction,
                                    forward, backward)
            enzymes.enzymes.append(enzyme)
            # machinery composition
            reactants = enzyme.machinery_composition.reactants
            for gene in comp:
                ref = self.data.protein_reference.get(gene, None)
                if ref:
                    reactants.append(rba.xml.SpeciesReference(*ref))
        # atpm enzyme
        reaction = self.default.atpm_reaction
        forward = backward = def_act.efficiency_id
        enzyme = rba.xml.Enzyme(reaction + '_enzyme', reaction,
                                forward, backward)
        enzymes.enzymes.append(enzyme)
        return enzymes

    def build_processes(self):
        """
        Build process part of RBA model.

        Returns
        -------
        rba.xml.RbaProcesses
            RBA process model in XML format.

        """
        processes = rba.xml.RbaProcesses()
        def_proc = DefaultProcesses(self.default, self.data.metabolite_map)
        compartments = self.data.compartments()
        # processes
        proc_list = processes.processes
        proc_list.append(def_proc.translation(
            {m.id: m.stoichiometry for m in self.data.ribosome}, compartments
            ))
        proc_list.append(def_proc.folding(
            {m.id: m.stoichiometry for m in self.data.chaperone}
            ))
        proc_list.append(def_proc.transcription())
        proc_list.append(def_proc.replication())
        proc_list.append(def_proc.rna_degradation())
        proc_list.append(def_proc.metabolite_production())
        proc_list.append(def_proc.macrocomponents(self.data.macrocomponents))
        proc_list.append(def_proc.maintenance_atp(self.default.atpm_reaction))
        # component maps
        map_list = processes.component_maps
        map_list.append(def_proc.translation_map(self.data.cofactors()))
        map_list.append(def_proc.folding_map())
        map_list.append(def_proc.transcription_map())
        map_list.append(def_proc.rna_degradation_map())
        map_list.append(def_proc.replication_map())
        return processes
