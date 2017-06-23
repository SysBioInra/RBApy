"""
Module defining DefaultProcesses class.
"""

# python 2/3 compatibility
from __future__ import division, print_function

# local imports
from rba import rba_xml

class DefaultProcesses(object):
    """
    Class initializing default process structure used by RBA.

    Attributes:
        rba_process: RBA process structure containing default processes.
        aggregates: xml structure containing a list of aggregates used by the
            default processes.
    """

    def __init__(self, default_data, ribosome_composition,
                 chaperone_composition, compartments, cofactors,
                 metabolite_map, macro_flux):
        """
        Constructor.

        Args:
            default_data (default_data.DefaultData): object containing default
                data used to initialize processes.
            ribosome_composition: dict mapping protein or rna identifiers to
                their stoichiometry within a ribosome unit.
            chaperone_composition: dict mapping protein or rna identifiers to
                their stoichiometry within an average chaperone.
            compartments: list of compartments in the model.
            cofactors (list of cofactors.Cofactor): cofactors in the model.
            metabolite_map: dict linking internal names of metabolites with
                user-defined SBML identifiers.
            macro_flux: dict mapping metabolite names with a production flux.
        """
        self._cofactors = cofactors
        self._metabolites = metabolite_map
        self._macro_flux = macro_flux
        self._target_already_treated = []

        self.rba_processes = rba_xml.RbaProcesses()
        self.aggregates = rba_xml.ListOfAggregates()

        # create default processes
        default_metabolites = default_data.metabolites
        default_parameters = default_data.parameters
        processes = self.rba_processes.processes
        processes.append(self._translation(default_metabolites,
                                           default_parameters,
                                           ribosome_composition,
                                           compartments))
        processes.append(self._folding(chaperone_composition))
        processes.append(self._transcription(default_metabolites))
        processes.append(self._replication(default_metabolites))
        processes.append(self._rna_degradation(default_metabolites))
        processes.append(self._metabolite_production())
        processes.append(self._macrocomponents())
        processes.append(self._maintenance_atp(default_data.atpm_reaction))

        # create default component maps
        component_maps = self.rba_processes.component_maps
        component_maps.append(self._translation_map(default_metabolites))
        component_maps.append(self._folding_map(default_metabolites))
        component_maps.append(self._transcription_map(default_metabolites))
        component_maps.append(self._rna_degradation_map(default_metabolites))
        component_maps.append(self._replication_map(default_metabolites))

    def _translation(self, default_metabolites, default_parameters,
                     ribosome_composition, compartments):
        process = rba_xml.Process('P_TA', 'Translation')
        # machinery
        machine = process.machinery.machinery_composition
        for id_, sto in ribosome_composition.items():
            machine.reactants.append(rba_xml.SpeciesReference(id_, sto))
        self._append_metabolite(machine.reactants, 'GTP', 2)
        self._append_metabolite(machine.reactants, 'H2O', 2)
        self._append_metabolite(machine.products, 'GDP', 2)
        self._append_metabolite(machine.products, 'Pi', 2)
        self._append_metabolite(machine.products, 'H', 2)
        id_ = 'ribosome_capacity'
        self._add_aggregate(id_, ['ribosome_efficiency_MM',
                                  'fraction_active_ribosomes'])
        process.machinery.capacity.value = id_
        # operating costs
        operation = rba_xml.Operation('translation', 'protein')
        process.operations.productions.append(operation)
        # targets
        for cpt in compartments:
            prot_id = default_metabolites.average_protein_id(cpt)
            fraction_fn = default_parameters.protein_fraction_id(cpt)
            non_enzymatic_fn = default_parameters.non_enzymatic_fraction_id(cpt)
            target = rba_xml.TargetSpecies(prot_id)
            id_ = 'nonenzymatic_proteins_' + cpt
            self._add_aggregate(id_, ['amino_acid_concentration',
                                      'inverse_average_protein_length',
                                      fraction_fn, non_enzymatic_fn])
            target.value = id_
            process.targets.concentrations.append(target)
        return process

    @staticmethod
    def _folding(chaperone_composition):
        process = rba_xml.Process('P_CHP', 'Folding')
        # capacity constraint
        machine = process.machinery.machinery_composition
        for id_, sto in chaperone_composition.items():
            machine.reactants.append(rba_xml.SpeciesReference(id_, sto))
        process.machinery.capacity.value = 'chaperone_efficiency_LM'
        # operating costs
        operation = rba_xml.Operation('folding', 'protein')
        process.operations.productions.append(operation)
        return process

    def _transcription(self, default_metabolites):
        process = rba_xml.Process('P_TSC', 'Transcription')
        # operating costs
        operation = rba_xml.Operation('transcription', 'rna')
        process.operations.productions.append(operation)
        # targets
        # mrna
        target = rba_xml.TargetSpecies(default_metabolites.mrna)
        target.value = 0.01
        process.targets.concentrations.append(target)
        target = rba_xml.TargetSpecies(default_metabolites.mrna)
        target.value = 0.15996
        process.targets.production_fluxes.append(target)
        # trnas
        for aa in default_metabolites.aas:
            key = default_metabolites.uncharged_trna_key(aa)
            self._target_already_treated.append(key)
            metabolite = self._metabolites[key]
            if metabolite.sbml_id and metabolite.concentration:
                target = rba_xml.TargetSpecies(metabolite.sbml_id)
                target.value = metabolite.concentration
                process.targets.concentrations.append(target)
        return process

    @staticmethod
    def _replication(default_metabolites):
        process = rba_xml.Process('P_REP', 'Replication')
        # operating costs
        operation = rba_xml.Operation('replication', 'dna')
        process.operations.productions.append(operation)
        # targets
        target = rba_xml.TargetSpecies(default_metabolites.dna)
        target.value = 0.0807
        process.targets.concentrations.append(target)
        return process

    @staticmethod
    def _rna_degradation(default_metabolites):
        process = rba_xml.Process('P_RNADEG', 'RNA degradation')
        # operating costs
        operation = rba_xml.Operation('rna_degradation', 'rna')
        process.operations.degradations.append(operation)
        # targets
        target = rba_xml.TargetSpecies(default_metabolites.mrna)
        target.value = 0.15996
        process.targets.degradation_fluxes.append(target)
        return process

    def _metabolite_production(self):
        process = rba_xml.Process('P_MET_PROD', 'Metabolite production')
        # targets
        for key in self._metabolites:
            if key in self._target_already_treated:
                continue
            metabolite = self._metabolites[key]
            # if a metabolite could not be identified, ignore it
            if metabolite.sbml_id and metabolite.concentration:
                target = rba_xml.TargetSpecies(metabolite.sbml_id)
                target.value = metabolite.concentration
                process.targets.concentrations.append(target)
        return process

    def _macrocomponents(self):
        process = rba_xml.Process('P_MACRO_PROD', 'Macrocomponent production')
        for id_, flux in self._macro_flux.items():
            target = rba_xml.TargetSpecies(id_)
            target.value = flux
            process.targets.concentrations.append(target)
        return process

    @staticmethod
    def _maintenance_atp(reaction_name):
        process = rba_xml.Process('P_maintenance_atp', 'Maintenance ATP')
        target = rba_xml.TargetReaction(reaction_name)
        target.lower_bound = 'maintenance_atp'
        process.targets.reaction_fluxes.append(target)
        return process

    def _translation_map(self, default_metabolites):
        map_ = rba_xml.ComponentMap('translation')

        # constant costs
        reactants = map_.constant_cost.reactants
        self._append_metabolite(reactants,
                                default_metabolites.charged_trna_key('fM'), 1)
        self._append_metabolite(reactants, 'GTP', 1)
        self._append_metabolite(reactants, 'H2O', 2)
        products = map_.constant_cost.products
        self._append_metabolite(products,
                                default_metabolites.uncharged_trna_key('M'), 1)
        self._append_metabolite(products, 'MET', 1)
        self._append_metabolite(products, 'FOR', 1)
        self._append_metabolite(products, 'GDP', 1)
        self._append_metabolite(products, 'Pi', 1)
        self._append_metabolite(products, 'H', 1)

        # amino acids
        for aa in default_metabolites.aas:
            cost = rba_xml.Cost(aa, 1)
            self._append_metabolite(cost.reactants,
                                    default_metabolites.charged_trna_key(aa), 1)
            self._append_metabolite(cost.reactants, 'GTP', 2)
            self._append_metabolite(cost.reactants, 'H2O', 2)
            self._append_metabolite(cost.products,
                                    default_metabolites.uncharged_trna_key(aa),
                                    1)
            self._append_metabolite(cost.products, 'GDP', 2)
            self._append_metabolite(cost.products, 'Pi', 2)
            self._append_metabolite(cost.products, 'H', 2)
            map_.costs.append(cost)

        # cofactors
        for cofactor in self._cofactors:
            cost = rba_xml.Cost(cofactor.id)
            self._append_metabolite(cost.reactants, cofactor.id, 1)
            map_.costs.append(cost)
        return map_

    @staticmethod
    def _folding_map(default_metabolites):
        map_ = rba_xml.ComponentMap('folding')
        for aa in default_metabolites.aas:
            map_.costs.append(rba_xml.Cost(aa, 0.1))
        return map_

    def _transcription_map(self, default_metabolites):
        map_ = rba_xml.ComponentMap('transcription')
        for n in default_metabolites.nucleotides:
            cost = rba_xml.Cost(n)
            self._append_metabolite(cost.reactants,
                                    default_metabolites.ntp_key(n), 1)
            self._append_metabolite(cost.reactants, 'H2O', 1)
            self._append_metabolite(cost.products, 'PPi', 1)
            self._append_metabolite(cost.products, 'H', 1)
            map_.costs.append(cost)
        return map_

    def _rna_degradation_map(self, default_metabolites):
        map_ = rba_xml.ComponentMap('rna_degradation')
        for n in default_metabolites.nucleotides:
            cost = rba_xml.Cost(n)
            self._append_metabolite(cost.reactants, 'H2O', 1)
            self._append_metabolite(cost.products,
                                    default_metabolites.nmp_key(n), 1)
            self._append_metabolite(cost.products, 'H', 1)
            map_.costs.append(cost)
        return map_

    def _replication_map(self, default_metabolites):
        map_ = rba_xml.ComponentMap('replication')
        for n in default_metabolites.d_nucleotides:
            cost = rba_xml.Cost(n)
            self._append_metabolite(cost.reactants,
                                    default_metabolites.dntp_key(n), 1)
            self._append_metabolite(cost.reactants, 'H2O', 1)
            self._append_metabolite(cost.products, 'PPi', 1)
            self._append_metabolite(cost.products, 'H', 1)
            map_.costs.append(cost)
        return map_

    def _append_metabolite(self, sr_list, key, sto):
        """
        Append species reference only if it was mapped to an SBML id.
        """
        sbml_id = self._metabolites[key].sbml_id
        if sbml_id:
            sr_list.append(rba_xml.SpeciesReference(sbml_id, sto))

    def _add_aggregate(self, id_, fn_refs):
        result = rba_xml.Aggregate(id_, 'multiplication')
        for ref in fn_refs:
            result.function_references.append(rba_xml.FunctionReference(ref))
        self.aggregates.append(result)

