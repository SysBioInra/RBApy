"""
Module defining DefaultProcesses class.
"""

# python 2/3 compatibility
from __future__ import division, print_function

# local imports
import rba.xml


class DefaultProcesses(object):
    """
    Class initializing default process structure used by RBA.

    Attributes:
        rba_process: RBA process structure containing default processes.
        aggregates: xml structure containing a list of aggregates used by the
            default processes.
    """

    def __init__(self, default_data, metabolite_map):
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
        self.default = default_data
        self._metabolites = metabolite_map

    def translation(self, ribosome_composition, compartments):
        process = rba.xml.Process('P_TA', 'Translation')
        # machinery
        for id_, sto in ribosome_composition.items():
            machine = process.machinery.machinery_composition
            machine.reactants.append(rba.xml.SpeciesReference(id_, sto))
        self._append_metabolite(machine.reactants, 'GTP', 2)
        self._append_metabolite(machine.reactants, 'H2O', 2)
        self._append_metabolite(machine.products, 'GDP', 2)
        self._append_metabolite(machine.products, 'Pi', 2)
        self._append_metabolite(machine.products, 'H', 2)
        process.machinery.capacity.value = 'ribosome_capacity'
        # operating costs
        operation = rba.xml.Operation('translation', 'protein')
        process.operations.productions.append(operation)
        # targets
        for cpt in compartments:
            prot_id = self.default.metabolites.average_protein_id(cpt)
            target = rba.xml.TargetSpecies(prot_id)
            target.value = 'nonenzymatic_proteins_' + cpt
            process.targets.concentrations.append(target)
        return process

    @staticmethod
    def folding(chaperone_composition):
        process = rba.xml.Process('P_CHP', 'Folding')
        # capacity constraint
        machine = process.machinery.machinery_composition
        for id_, sto in chaperone_composition.items():
            machine.reactants.append(rba.xml.SpeciesReference(id_, sto))
        process.machinery.capacity.value = 'chaperone_efficiency_LM'
        # operating costs
        operation = rba.xml.Operation('folding', 'protein')
        process.operations.productions.append(operation)
        return process

    def transcription(self):
        process = rba.xml.Process('P_TSC', 'Transcription')
        default_metabolites = self.default.metabolites
        # operating costs
        operation = rba.xml.Operation('transcription', 'rna')
        process.operations.productions.append(operation)
        # targets
        # mrna
        target = rba.xml.TargetSpecies(default_metabolites.mrna)
        target.value = 0.01
        process.targets.concentrations.append(target)
        target = rba.xml.TargetSpecies(default_metabolites.mrna)
        target.value = 0.15996
        process.targets.production_fluxes.append(target)
        # trnas
        for aa in default_metabolites.aas:
            key = default_metabolites.uncharged_trna_key(aa)
            metabolite = self._metabolites[key]
            if metabolite.sbml_id and metabolite.concentration:
                target = rba.xml.TargetSpecies(metabolite.sbml_id)
                target.value = metabolite.concentration
                process.targets.concentrations.append(target)
        return process

    def replication(self):
        process = rba.xml.Process('P_REP', 'Replication')
        # operating costs
        operation = rba.xml.Operation('replication', 'dna')
        process.operations.productions.append(operation)
        # targets
        target = rba.xml.TargetSpecies(self.default.metabolites.dna)
        target.value = 0.0807
        process.targets.concentrations.append(target)
        return process

    def rna_degradation(self):
        process = rba.xml.Process('P_RNADEG', 'RNA degradation')
        # operating costs
        operation = rba.xml.Operation('rna_degradation', 'rna')
        process.operations.degradations.append(operation)
        # targets
        target = rba.xml.TargetSpecies(self.default.metabolites.mrna)
        target.value = 0.15996
        process.targets.degradation_fluxes.append(target)
        return process

    def metabolite_production(self):
        process = rba.xml.Process('P_MET_PROD', 'Metabolite production')
        # targets
        def_metab = self.default.metabolites
        uncharged_trnas = [def_metab.uncharged_trna_key(aa)
                           for aa in def_metab.aas]
        for key, metabolite in self._metabolites.items():
            if key in uncharged_trnas:
                continue
            # if a metabolite could not be identified, ignore it
            if metabolite.sbml_id and metabolite.concentration:
                target = rba.xml.TargetSpecies(metabolite.sbml_id)
                target.value = metabolite.concentration
                process.targets.concentrations.append(target)
        return process

    @staticmethod
    def macrocomponents(macro_fluxes):
        process = rba.xml.Process('P_MACRO_PROD', 'Macrocomponent production')
        for id_, flux in macro_fluxes.items():
            target = rba.xml.TargetSpecies(id_)
            target.value = flux
            process.targets.concentrations.append(target)
        return process

    @staticmethod
    def maintenance_atp(reaction_name):
        process = rba.xml.Process('P_maintenance_atp', 'Maintenance ATP')
        target = rba.xml.TargetReaction(reaction_name)
        target.lower_bound = 'maintenance_atp'
        process.targets.reaction_fluxes.append(target)
        return process

    def translation_map(self, cofactors):
        map_ = rba.xml.ComponentMap('translation')
        def_metabolites = self.default.metabolites
        # constant costs
        reactants = map_.constant_cost.reactants
        self._append_metabolite(
            reactants, def_metabolites.charged_trna_key('fM'), 1
            )
        self._append_metabolite(reactants, 'GTP', 1)
        self._append_metabolite(reactants, 'H2O', 2)
        products = map_.constant_cost.products
        self._append_metabolite(products,
                                def_metabolites.uncharged_trna_key('M'), 1)
        self._append_metabolite(products, 'MET', 1)
        self._append_metabolite(products, 'FOR', 1)
        self._append_metabolite(products, 'GDP', 1)
        self._append_metabolite(products, 'Pi', 1)
        self._append_metabolite(products, 'H', 1)
        # amino acids
        for aa in def_metabolites.aas:
            cost = rba.xml.Cost(aa, 1)
            self._append_metabolite(
                cost.reactants, def_metabolites.charged_trna_key(aa), 1
                )
            self._append_metabolite(cost.reactants, 'GTP', 2)
            self._append_metabolite(cost.reactants, 'H2O', 2)
            self._append_metabolite(
                cost.products, def_metabolites.uncharged_trna_key(aa), 1
                )
            self._append_metabolite(cost.products, 'GDP', 2)
            self._append_metabolite(cost.products, 'Pi', 2)
            self._append_metabolite(cost.products, 'H', 2)
            map_.costs.append(cost)
        # cofactors
        for cofactor in cofactors:
            cost = rba.xml.Cost(cofactor.chebi)
            self._append_metabolite(cost.reactants, cofactor.chebi, 1)
            map_.costs.append(cost)
        return map_

    def folding_map(self):
        map_ = rba.xml.ComponentMap('folding')
        for aa in self.default.metabolites.aas:
            map_.costs.append(rba.xml.Cost(aa, 0.1))
        return map_

    def transcription_map(self):
        map_ = rba.xml.ComponentMap('transcription')
        for n in self.default.metabolites.nucleotides:
            cost = rba.xml.Cost(n)
            self._append_metabolite(
                cost.reactants, self.default.metabolites.ntp_key(n), 1
                )
            self._append_metabolite(cost.reactants, 'H2O', 1)
            self._append_metabolite(cost.products, 'PPi', 1)
            self._append_metabolite(cost.products, 'H', 1)
            map_.costs.append(cost)
        return map_

    def rna_degradation_map(self):
        map_ = rba.xml.ComponentMap('rna_degradation')
        for n in self.default.metabolites.nucleotides:
            cost = rba.xml.Cost(n)
            self._append_metabolite(cost.reactants, 'H2O', 1)
            self._append_metabolite(
                cost.products, self.default.metabolites.nmp_key(n), 1
                )
            self._append_metabolite(cost.products, 'H', 1)
            map_.costs.append(cost)
        return map_

    def replication_map(self):
        map_ = rba.xml.ComponentMap('replication')
        for n in self.default.metabolites.d_nucleotides:
            cost = rba.xml.Cost(n)
            self._append_metabolite(
                cost.reactants, self.default.metabolites.dntp_key(n), 1
                )
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
            sr_list.append(rba.xml.SpeciesReference(sbml_id, sto))


def aggregates(default_params, compartments):
    result = []
    result.append(build_aggregate(
        'ribosome_capacity',
        ['ribosome_efficiency_MM', 'fraction_active_ribosomes']
        ))
    for cpt in compartments:
        fraction_fn = default_params.protein_fraction_id(cpt)
        non_enzymatic_fn = default_params.non_enzymatic_fraction_id(cpt)
        result.append(build_aggregate(
            'nonenzymatic_proteins_' + cpt,
            ['amino_acid_concentration', 'inverse_average_protein_length',
             fraction_fn, non_enzymatic_fn]
            ))
    return result


def build_aggregate(id_, fn_refs):
    result = rba.xml.Aggregate(id_, 'multiplication')
    for ref in fn_refs:
        result.function_references.append(rba.xml.FunctionReference(ref))
    return result
