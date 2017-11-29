"""Module defining DefaultProcesses class."""

# python 2/3 compatibility
from __future__ import division, print_function, absolute_import

# local imports
import rba.xml


class DefaultProcesses(object):
    """
    Class initializing default process structure used by RBA.

    Attributes
    ----------
    default : rba.prerba.default_data.DefaultData
        Default data.

    """

    def __init__(self, default_data, metabolite_map):
        """
        Build object from default data and metabolite map.

        Parameters
        ----------
        default_data : rba.prerba.default_data.DefaultData)
            Default data used to initialize processes.
        metabolite_map: dict
            Map from internal metabolite ids to user-defined ids.

        """
        self.default = default_data
        self._metabolites = metabolite_map

    def translation(self, ribosome_composition):
        """
        Build translation process.

        Parameters
        ----------
        ribosome_composition : dict
            Map from molecule identifiers to
            their stoichiometry within a ribosome unit.

        Returns
        -------
        rba.xml.Process
            Translation process.

        """
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
        return process

    def folding(self, chaperone_composition):
        """
        Build folding process.

        Parameters
        ----------
        chaperone_composition : dict
            Map from molecule identifiers to
            their stoichiometry within an average chaperone.

        Returns
        -------
        rba.xml.Process
            Folding process.

        """
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
        """
        Build transcription process.

        Returns
        -------
        rba.xml.Process
            Transcription process.

        """
        process = rba.xml.Process('P_TSC', 'Transcription')
        default_metabolites = self.default.metabolites
        # operating costs
        operation = rba.xml.Operation('transcription', 'rna')
        process.operations.productions.append(operation)
        return process

    def replication(self):
        """
        Build replication process.

        Returns
        -------
        rba.xml.Process
            Replication process.

        """
        process = rba.xml.Process('P_REP', 'Replication')
        # operating costs
        operation = rba.xml.Operation('replication', 'dna')
        process.operations.productions.append(operation)
        return process

    def rna_degradation(self):
        """
        Build RNA degradation process.

        Returns
        -------
        rba.xml.Process
            RNA degradation process.

        """
        process = rba.xml.Process('P_RNADEG', 'RNA degradation')
        # operating costs
        operation = rba.xml.Operation('rna_degradation', 'rna')
        process.operations.degradations.append(operation)
        return process

    def translation_map(self, cofactors):
        """
        Build translation map.

        Parameters
        ----------
        cofactors : list of rba.prerba.uniprot_data.Cofactor
            Cofactors in the model.

        Returns
        -------
        rba.xml.ComponentMap
            Translation map.

        """
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
        """
        Build folding map.

        Returns
        -------
        rba.xml.ComponentMap
            Folding map.

        """
        map_ = rba.xml.ComponentMap('folding')
        for aa in self.default.metabolites.aas:
            map_.costs.append(rba.xml.Cost(aa, 0.1))
        return map_

    def transcription_map(self):
        """
        Build transcription map.

        Returns
        -------
        rba.xml.ComponentMap
            Transcription map.

        """
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
        """
        Build RNA degradation map.

        Returns
        -------
        rba.xml.ComponentMap
            RNA degradation map.

        """
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
        """
        Build replication map.

        Returns
        -------
        rba.xml.ComponentMap
            Replication map.

        """
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
        """Append species reference only if it was mapped to an SBML id."""
        sbml_id = self._metabolites[key].sbml_id
        if sbml_id:
            sr_list.append(rba.xml.SpeciesReference(sbml_id, sto))
