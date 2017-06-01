
import xml.dom.minidom
import copy

from ..rba_xml import *
from . import rba_data

class DefaultProcesses(object):
    """
    """
    def __init__(self, ribosome_composition, chaperone_composition,
                 compartments, cofactors, metabolites, macro_flux):
        self._cofactors = cofactors
        self._metabolites = metabolites
        self._macro_flux = macro_flux
        self._target_already_treated = []
        self._missing_id = ''
        
        self.rba_processes = RbaProcesses()
        self.aggregates = ListOfAggregates()

        # create default processes
        processes = self.rba_processes.processes
        processes.append(self._translation(ribosome_composition, compartments))
        processes.append(self._folding(chaperone_composition))
        processes.append(self._transcription())
        processes.append(self._replication())
        processes.append(self._rna_degradation())
        processes.append(self._metabolite_production())
        processes.append(self._macrocomponents())
        processes.append(self._maintenance_atp())

        # create default component maps
        component_maps = self.rba_processes.component_maps
        component_maps.append(self._translation_map())
        component_maps.append(self._folding_map())
        component_maps.append(self._transcription_map())
        component_maps.append(self._rna_degradation_map())
        component_maps.append(self._replication_map())

    def _translation(self, ribosome_composition, compartments):
        process = Process('P_TA', 'Translation')
        # machinery
        machine = process.machinery.machinery_composition
        for id_, sto in ribosome_composition.items():
            machine.reactants.append(SpeciesReference(id_, sto))
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
        process.operations.productions\
                               .append(Operation('translation', 'protein'))
        # targets
        for c in compartments:
            target = TargetSpecies(rba_data.average_protein_id(c))
            id_ = 'nonenzymatic_proteins_' + c
            self._add_aggregate(id_, ['amino_acid_concentration',
                                      'inverse_average_protein_length',
                                      rba_data.protein_fraction_id(c),
                                      rba_data.non_enzymatic_fraction_id(c)])
            target.value = id_
            process.targets.concentrations.append(target)
        return process

    def _folding(self, chaperone_composition):
        process = Process('P_CHP', 'Folding')
        # capacity constraint
        machine = process.machinery.machinery_composition
        for id_, sto in chaperone_composition.items():
            machine.reactants.append(SpeciesReference(id_, sto))
        process.machinery.capacity.value = 'chaperone_efficiency_LM'
        # operating costs
        process.operations.productions\
                               .append(Operation('folding', 'protein'))
        return process

    def _transcription(self):
        process = Process('P_TSC', 'Transcription')
        # operating costs
        process.operations.productions\
                               .append(Operation('transcription', 'rna'))
        # targets
        # mrna
        target = TargetSpecies(rba_data.mrna)
        target.value = 0.01
        process.targets.concentrations.append(target)
        target = TargetSpecies(rba_data.mrna)
        target.value = 0.15996
        process.targets.production_fluxes.append(target)
        # trnas
        for aa in rba_data.aas:
            key = rba_data.uncharged_trna_key(aa)
            self._target_already_treated.append(key)
            id_ = self._metabolites[key].sbml_id
            if id_ == self._missing_id: continue
            target = TargetSpecies(id_)
            target.value = self._metabolites[key].concentration
            process.targets.concentrations.append(target)
        return process

    def _replication(self):
        process = Process('P_REP', 'Replication')
        # operating costs
        process.operations.productions\
                               .append(Operation('replication', 'dna'))
        # targets
        target = TargetSpecies(rba_data.dna)
        target.value = 0.0807
        process.targets.concentrations.append(target)
        return process

    def _rna_degradation(self):
        process = Process('P_RNADEG', 'RNA degradation')
        # operating costs
        process.operations.degradations\
                               .append(Operation('rna_degradation', 'rna'))
        # targets
        target = TargetSpecies(rba_data.mrna)
        target.value = 0.15996
        process.targets.degradation_fluxes.append(target)
        return process

    def _metabolite_production(self):
        process = Process('P_MET_PROD', 'Metabolite production')
        # targets
        for key in self._metabolites:
            if key in self._target_already_treated: continue
            id_ = self._metabolites[key].sbml_id
            # if a metabolite could not be identified, ignore it
            if id_ == self._missing_id: continue
            flux = self._metabolites[key].concentration
            if flux == 0: continue
            target = TargetSpecies(id_)
            target.value = flux
            process.targets.concentrations.append(target)
        return process

    def _macrocomponents(self):
        process = Process('P_MACRO_PROD', 'Macrocomponent production')
        for id_, flux in self._macro_flux.items():
            target = TargetSpecies(id_)
            target.value = flux
            process.targets.concentrations.append(target)
        return process
    
    def _maintenance_atp(self):
        process = Process('P_maintenance_atp', 'Maintenance ATP')
        target = TargetReaction(rba_data.atpm_reaction)
        target.lower_bound = 'maintenance_atp'
        process.targets.reaction_fluxes.append(target)
        return process

    def _translation_map(self):
        map_ = ComponentMap('translation')
        
        # constant costs
        reactants = map_.constant_cost.reactants
        self._append_metabolite(reactants, rba_data.charged_trna_key('fM'), 1)
        self._append_metabolite(reactants, 'GTP', 1)
        self._append_metabolite(reactants, 'H2O', 2)
        products = map_.constant_cost.products
        self._append_metabolite(products, rba_data.uncharged_trna_key('M'), 1)
        self._append_metabolite(products, 'MET', 1)
        self._append_metabolite(products, 'FOR', 1)
        self._append_metabolite(products, 'GDP', 1)
        self._append_metabolite(products, 'Pi', 1)        
        self._append_metabolite(products, 'H', 1)

        # amino acids
        for aa in rba_data.aas:
            cost = Cost(aa, 1)
            self._append_metabolite(cost.reactants,
                                    rba_data.charged_trna_key(aa), 1)
            self._append_metabolite(cost.reactants, 'GTP', 2)
            self._append_metabolite(cost.reactants, 'H2O', 2)
            self._append_metabolite(cost.products,
                                    rba_data.uncharged_trna_key(aa), 1)
            self._append_metabolite(cost.products, 'GDP', 2)
            self._append_metabolite(cost.products, 'Pi', 2)        
            self._append_metabolite(cost.products, 'H', 2)
            map_.costs.append(cost)

        # cofactors
        for c in self._cofactors:
            cost = Cost(c.id)
            self._append_metabolite(cost.reactants, c.id, 1)
            map_.costs.append(cost)
        return map_

    def _folding_map(self):
        map_ = ComponentMap('folding')        
        for aa in rba_data.aas:
            map_.costs.append(Cost(aa, 0.1))
        return map_

    def _transcription_map(self):
        map_ = ComponentMap('transcription')        
        for n in rba_data.nucleotides:
            cost = Cost(n)
            self._append_metabolite(cost.reactants, rba_data.ntp_key(n), 1)
            self._append_metabolite(cost.reactants, 'H2O', 1)
            self._append_metabolite(cost.products, 'PPi', 1)
            self._append_metabolite(cost.products, 'H', 1)
            map_.costs.append(cost)
        return map_

    def _rna_degradation_map(self):
        map_ = ComponentMap('rna_degradation')
        for n in rba_data.nucleotides:
            cost = Cost(n)
            self._append_metabolite(cost.reactants, 'H2O', 1)
            self._append_metabolite(cost.products, rba_data.nmp_key(n), 1)
            self._append_metabolite(cost.products, 'H', 1)
            map_.costs.append(cost)
        return map_

    def _replication_map(self):
        map_ = ComponentMap('replication')
        for n in rba_data.d_nucleotides:
            cost = Cost(n)
            self._append_metabolite(cost.reactants, rba_data.dntp_key(n), 1)
            self._append_metabolite(cost.reactants, 'H2O', 1)
            self._append_metabolite(cost.products, 'PPi', 1)
            self._append_metabolite(cost.products, 'H', 1)
            map_.costs.append(cost)
        return map_        

    def _append_metabolite(self, SR_list, key, sto):
        # append species reference only if a corresponding SBML id has
        # been defined.
        sbml_id = self._metabolites[key].sbml_id
        if sbml_id != self._missing_id:
            SR_list.append(SpeciesReference(sbml_id, sto))

    def _add_aggregate(self, id_, fn_refs):
        result = Aggregate(id_, 'multiplication')
        for ref in fn_refs:
            result.function_references.append(FunctionReference(ref))
        self.aggregates.append(result)

