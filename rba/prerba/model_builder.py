"""Model used to build RBA XML model."""

# python 2/3 compatibility
from __future__ import division, print_function, absolute_import

# global imports
import itertools
import copy

# local imports
from rba import RbaModel
from rba.prerba.user_data import UserData
from rba.prerba.default_data import DefaultData
from rba.prerba.user_data import ntp_composition
from rba.prerba.default_processes import DefaultProcesses
from rba.prerba.default_targets import DefaultTargets
import rba.xml


class ModelBuilder(object):
    """Build a RBA model from user data."""

    def __init__(self, parameter_file):
        """Constructor."""
        self.data = UserData(parameter_file)
        self.default = DefaultData()

    def build_model(self):
        """Build and return entire RbaModel."""
        model = RbaModel()
        model.metabolism = self.build_metabolism()
        model.density = self.build_density()
        model.parameters = self.build_parameters()
        model.proteins = self.build_proteins()
        model.rnas = self.build_rnas()
        model.dna = self.build_dna()
        model.processes = self.build_processes()
        model.targets = self.build_targets()
        model.enzymes = self.build_enzymes()
        model.medium = self.build_medium()
        model.output_dir = self.data.output_dir()
        return model

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
        metabolism.reactions.append(self._atpm_reaction())
        for c in self.data.compartments():
            metabolism.compartments.append(rba.xml.Compartment(c))
        return metabolism

    def _atpm_reaction(self):
        reaction = rba.xml.Reaction(self.default.atpm_reaction, False)
        for m in ['ATP', 'H2O']:
            id_ = self.data.metabolite_map[m].sbml_id
            if id_:
                reaction.reactants.append(rba.xml.SpeciesReference(id_, 1))
        for m in ['ADP', 'H', 'Pi']:
            id_ = self.data.metabolite_map[m].sbml_id
            if id_:
                reaction.products.append(rba.xml.SpeciesReference(id_, 1))
        return reaction

    def build_density(self):
        """
        Build density part of RBA model.

        Returns
        -------
        rba.xml.RbaDensity
            RBA density model in XML format.

        """
        density = rba.xml.RbaDensity()
        external = self._external()
        for c_id in self.data.compartments():
            if c_id != external:
                new_density = rba.xml.TargetDensity(c_id)
                new_density.upper_bound = c_id + '_density'
                density.target_densities.append(new_density)
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
        for fn in self._all_parameter_functions():
            parameters.functions.append(fn)
        for fn in self._all_parameter_aggregates():
            parameters.aggregates.append(agg)
        return parameters

    def _all_parameter_functions(self):
        return (self._density_functions() + self._protein_functions()
                + self.default.parameters.process_functions()
                + self._target_functions() + self._efficiency_functions())

    def _density_functions(self):
        return self.default.parameters.density_functions(
            self._cytoplasm(), self._external(), self._other_compartments()
            )

    def _cytoplasm(self):
        return self.data.compartment('Cytoplasm')

    def _external(self):
        return self.data.compartment('Secreted')

    def _other_compartments(self):
        result = self.data.compartments()
        result.remove(self._cytoplasm())
        result.remove(self._external())
        return result

    def _protein_functions(self):
        return [self.default.parameters.inverse_average_protein_length(
            self.data.average_protein_length()
        )]

    def _target_functions(self):
        return [
            self.default.parameters.metabolite_concentration_function(*target)
            for target in self.data.metabolite_targets()
        ]

    def _efficiency_functions(self):
        fns = [self.default.activity.efficiency_function(),
               self.default.activity.transport_function()]
        for r_id in self.data.transport_reaction_ids():
            fns += self.default.activity.transport_functions(
                r_id, self.data.imported_metabolites(r_id)
            )
        return fns

    def _append_all_parameter_aggregates(self, parameters):
        self._append_aggregates(parameters, self._density_aggregates())
        self._append_aggregates(parameters, self._process_aggregates())
        self._append_aggregates(parameters, self._efficiency_aggregates())

    def _append_aggregates(self, parameters, aggs):
        for agg in aggs:
            parameters.aggregates.append(agg)

    def _density_aggregates(self):
        return self.default.parameters.density_aggregates(
            self._cytoplasm(), self._external(), self._other_compartments()
            )

    def _process_aggregates(self):
        return self.default.parameters.process_aggregates()

    def _efficiency_aggregates(self):
        return [self.default.activity.transport_aggregate(
                    r, self.data.imported_metabolites(r)
                )
                for r in self.data.transport_reaction_ids()]

    def build_proteins(self):
        """
        Build protein part of RBA model.

        Returns
        -------
        rba.xml.RbaMacromolecules
            RBA protein model in XML format.

        """
        builder = MacromoleculeBuilder()
        # components
        for aa in self.default.metabolites.aas:
            builder.add_component(aa, '', 'amino_acid', 1)
        for c in self.data.cofactors():
            builder.add_component(c.chebi, c.name, 'cofactor', 0)
        # enzymatic proteins
        for gene_name, protein in self.data.enzymatic_proteins.items():
            builder.add_macromolecule(gene_name, protein.location,
                                      self._protein_composition(protein))
        # average proteins
        for comp in self.data.compartments():
            builder.add_macromolecule(self.data.average_protein_id(comp),
                                      comp, self.data.average_protein())
        # machinery proteins
        for prot in itertools.chain(self.data.ribosome.proteins,
                                    self.data.chaperone.proteins):
            builder.add_macromolecule(prot.id, self._cytoplasm(),
                                      self.data.aa_composition(prot.sequence))
        return builder.result

    def _protein_composition(self, protein):
        comp = self.data.aa_composition(protein.sequence)
        for cofactor in protein.cofactors:
            comp[cofactor.chebi] = cofactor.stoichiometry
        return comp

    def build_rnas(self):
        """
        Build RNA part of RBA model.

        Returns
        -------
        rba.xml.RbaMacromolecules
            RBA RNA model in XML format.

        """
        builder = MacromoleculeBuilder()
        builder.add_component('A', 'Adenosine residue', 'Nucleotide', 2.9036)
        builder.add_component('C', 'Cytosine residue', 'Nucleotide', 2.7017)
        builder.add_component('G', 'Guanine residue', 'Nucleotide', 3.0382)
        builder.add_component('U', 'Uramine residue', 'Nucleotide', 2.7102)
        # user rnas
        for rna_id, composition in self.data.rna_data.items():
            builder.add_macromolecule(rna_id, self._cytoplasm(), composition)
        # average RNA
        builder.add_macromolecule(
            self.default.metabolites.mrna, self._cytoplasm(),
            {'A': 0.2818, 'C': 0.2181, 'G': 0.2171, 'U': 0.283}
            )
        # machinery rnas
        for rna in itertools.chain(self.data.ribosome.rnas,
                                   self.data.chaperone.rnas):
            builder.add_macromolecule(rna.id, self._cytoplasm(),
                                      ntp_composition(rna.sequence))
        return builder.result

    def build_dna(self):
        """
        Build DNA part of RBA model.

        Returns
        -------
        rba.xml.RbaMacromolecules
            RBA DNA model in XML format.

        """
        builder = MacromoleculeBuilder()
        builder.add_component('A', 'Adenosine residue', 'Nucleotide', 0)
        builder.add_component('C', 'Cytosine residue', 'Nucleotide', 0)
        builder.add_component('G', 'Guanine residue', 'Nucleotide', 0)
        builder.add_component('T', 'Thymine residue', 'Nucleotide', 0)
        builder.add_macromolecule(
            self.default.metabolites.dna, self._cytoplasm(),
            {'A': 0.2818, 'C': 0.2181, 'G': 0.2171, 'T': 0.283}
            )
        return builder.result

    def build_enzymes(self):
        """
        Build enzyme part of RBA model.

        Returns
        -------
        rba.xml.RbaEnzymes
            RBA enzyme model in XML format.

        """
        enzymes = rba.xml.RbaEnzymes()
        reaction_ids = (r.id for r in self.data.sbml_reactions())
        for r_id, comp in zip(reaction_ids, self.data.enzyme_composition()):
            enzymes.enzymes.append(self._build_user_enzyme(r_id, comp))
        enzymes.enzymes.append(self._atpm_enzyme())
        return enzymes

    def _build_user_enzyme(self, reaction, composition):
        forward, backward = self._reaction_efficiencies(reaction)
        enzyme = rba.xml.Enzyme(reaction + '_enzyme', reaction,
                                forward, backward)
        self._add_enzyme_machinery(enzyme, composition)
        return enzyme

    def _reaction_efficiencies(self, reaction):
        if self.data.is_transport_reaction(reaction):
            forward = self.default.activity.transport_aggregate_id(reaction)
            backward = self.default.activity.transport_id
        else:
            forward = backward = self.default.activity.efficiency_id
        return forward, backward

    def _add_enzyme_machinery(self, enzyme, composition):
        # machinery composition
        reactants = enzyme.machinery_composition.reactants
        for gene in composition:
            ref = self.data.protein_reference.get(gene, None)
            if ref:
                reactants.append(rba.xml.SpeciesReference(*ref))

    def _atpm_enzyme(self):
        r_id = self.default.atpm_reaction
        forward = backward = self.default.activity.efficiency_id
        return rba.xml.Enzyme(r_id + '_enzyme', r_id, forward, backward)

    def build_processes(self):
        """
        Build process part of RBA model.

        Returns
        -------
        rba.xml.RbaProcesses
            RBA process model in XML format.

        """
        processes = rba.xml.RbaProcesses()
        for p in self._all_processes():
            processes.processes.append(p)
        for m in self._all_component_maps():
            processes.processing_maps.append(m)
        return processes

    def _all_processes(self):
        def_proc = DefaultProcesses(self.default, self.data.metabolite_map)
        proteins = self._all_protein_ids()
        rnas = self._all_rna_ids()
        dnas = self._all_dna_ids()
        result = [
            def_proc.translation(self.data.ribosome.composition(), proteins),
            def_proc.folding(self.data.chaperone.composition(), proteins),
            def_proc.transcription(rnas),
            def_proc.replication(dnas),
            def_proc.rna_degradation(rnas)
        ]
        # !!! Useless processes added because test model runs faster
        # with these empty processes
        result += [rba.xml.Process('test_process_{}'.format(i), 'Test process')
                   for i in range(3)]
        return result

    def _all_protein_ids(self):
        proteins = list(self.data.enzymatic_proteins.keys())
        proteins += [self.data.average_protein_id(c)
                     for c in self.data.compartments()]
        proteins += self.data.ribosome.protein_ids()
        proteins += self.data.chaperone.protein_ids()
        return proteins

    def _all_rna_ids(self):
        rnas = list(self.data.rna_data.keys())
        rnas.append(self.default.metabolites.mrna)
        rnas += self.data.ribosome.rna_ids()
        rnas += self.data.chaperone.rna_ids()
        return rnas

    def _all_dna_ids(self):
        return [self.default.metabolites.dna]

    def _all_component_maps(self):
        def_proc = DefaultProcesses(self.default, self.data.metabolite_map)
        return [
            def_proc.translation_map(self.data.cofactors()),
            def_proc.folding_map(),
            def_proc.transcription_map(),
            def_proc.rna_degradation_map(),
            def_proc.replication_map()
        ]

    def build_targets(self):
        """
        Build target part of RBA model.

        Returns
        -------
        rba.xml.RbaTargets
            RBA targets in XML format.

        """
        targets = rba.xml.RbaTargets()
        for t in self._all_targets():
            targets.target_groups.append(t)
        return targets

    def _all_targets(self):
        def_targ = DefaultTargets(self.default, self.data.metabolite_map)
        return [
                def_targ.translation(self.data.compartments()),
            def_targ.transcription(),
            def_targ.replication(),
            def_targ.rna_degradation(),
            def_targ.metabolite_production(),
            def_targ.macrocomponents(self.data.macrocomponents),
            def_targ.maintenance_atp(self.default.atpm_reaction)
        ]

    def build_medium(self):
        # !!! we identify metabolites by their prefix !!!
        # M_glc_p and M_glc_e will be seen as the same metabolite M_glc
        prefixes = (m.rsplit('_', 1)[0]
                    for m in self.data.external_metabolites())
        return dict.fromkeys(prefixes,
                             self.default.activity.medium_concentration)

    def export_proteins(self, filename):
        self.data.export_proteins(filename)


class MacromoleculeBuilder(object):
    def __init__(self):
        self.result = rba.xml.RbaMacromolecules()

    def add_component(self, id_, name, type_, stoichiometry):
        self.result.components.append(
            rba.xml.Component(id_, name, type_, stoichiometry)
            )

    def add_macromolecule(self, id_, location, composition):
        self.result.macromolecules.append(
            rba.xml.Macromolecule(id_, location, composition)
            )
