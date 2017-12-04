"""Module defining DefaultProcesses class."""

# python 2/3 compatibility
from __future__ import division, print_function, absolute_import

# local imports
import rba.xml


def create_processing(processing_map, set_, inputs):
    """Create rba.xml.Processing with given information."""
    processing = rba.xml.Processing(processing_map, set_)
    for id_ in inputs:
        processing.inputs.append(rba.xml.SpeciesReference(id_, 1))
    return processing


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

    def translation(self, ribosome_composition, inputs):
        """
        Build translation process.

        Parameters
        ----------
        ribosome_composition : dict
            Map from molecule identifiers to
            their stoichiometry within a ribosome unit.
        inputs : list of str
            Identifiers of process inputs.

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
        # processings
        process.processings.productions.append(
            create_processing('translation', 'protein', inputs)
            )
        return process

    def folding(self, chaperone_composition, inputs):
        """
        Build folding process.

        Parameters
        ----------
        chaperone_composition : dict
            Map from molecule identifiers to
            their stoichiometry within an average chaperone.
        inputs : list of str
            Identifiers of process inputs.

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
        process.processings.productions.append(
            create_processing('folding', 'protein', inputs)
            )
        return process

    def transcription(self, inputs):
        """
        Build transcription process.

        Parameters
        ----------
        inputs : list of str
            Identifiers of process inputs.

        Returns
        -------
        rba.xml.Process
            Transcription process.

        """
        process = rba.xml.Process('P_TSC', 'Transcription')
        default_metabolites = self.default.metabolites
        process.processings.productions.append(
            create_processing('transcription', 'rna', inputs)
            )
        return process

    def replication(self, inputs):
        """
        Build replication process.

        Parameters
        ----------
        inputs : list of str
            Identifiers of process inputs.

        Returns
        -------
        rba.xml.Process
            Replication process.

        """
        process = rba.xml.Process('P_REP', 'Replication')
        process.processings.productions.append(
            create_processing('replication', 'dna', inputs)
            )
        return process

    def rna_degradation(self, inputs):
        """
        Build RNA degradation process.

        Parameters
        ----------
        inputs : list of str
            Identifiers of process inputs.

        Returns
        -------
        rba.xml.Process
            RNA degradation process.

        """
        process = rba.xml.Process('P_RNADEG', 'RNA degradation')
        process.processings.degradations.append(
            create_processing('rna_degradation', 'rna', inputs)
            )
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
        rba.xml.ProcessingMap
            Translation map.

        """
        map_ = rba.xml.ProcessingMap('translation')
        def_metabolites = self.default.metabolites
        # constant costs
        reactants = map_.constant_processing.reactants
        self._append_metabolite(
            reactants, def_metabolites.charged_trna_key('fM'), 1
            )
        self._append_metabolite(reactants, 'GTP', 1)
        self._append_metabolite(reactants, 'H2O', 2)
        products = map_.constant_processing.products
        self._append_metabolite(products,
                                def_metabolites.uncharged_trna_key('M'), 1)
        self._append_metabolite(products, 'MET', 1)
        self._append_metabolite(products, 'FOR', 1)
        self._append_metabolite(products, 'GDP', 1)
        self._append_metabolite(products, 'Pi', 1)
        self._append_metabolite(products, 'H', 1)
        # amino acids
        for aa in def_metabolites.aas:
            cost = rba.xml.ComponentProcessing(aa, 1)
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
            map_.component_processings.append(cost)
        # cofactors
        for cofactor in cofactors:
            cost = rba.xml.ComponentProcessing(cofactor.chebi)
            self._append_metabolite(cost.reactants, cofactor.chebi, 1)
            map_.component_processings.append(cost)
        return map_

    def folding_map(self):
        """
        Build folding map.

        Returns
        -------
        rba.xml.ProcessingMap
            Folding map.

        """
        map_ = rba.xml.ProcessingMap('folding')
        for aa in self.default.metabolites.aas:
            map_.component_processings.append(
                rba.xml.ComponentProcessing(aa, 0.1)
                )
        return map_

    def transcription_map(self):
        """
        Build transcription map.

        Returns
        -------
        rba.xml.ProcessingMap
            Transcription map.

        """
        map_ = rba.xml.ProcessingMap('transcription')
        for n in self.default.metabolites.nucleotides:
            cost = rba.xml.ComponentProcessing(n)
            self._append_metabolite(
                cost.reactants, self.default.metabolites.ntp_key(n), 1
                )
            self._append_metabolite(cost.reactants, 'H2O', 1)
            self._append_metabolite(cost.products, 'PPi', 1)
            self._append_metabolite(cost.products, 'H', 1)
            map_.component_processings.append(cost)
        return map_

    def rna_degradation_map(self):
        """
        Build RNA degradation map.

        Returns
        -------
        rba.xml.ProcessingMap
            RNA degradation map.

        """
        map_ = rba.xml.ProcessingMap('rna_degradation')
        for n in self.default.metabolites.nucleotides:
            cost = rba.xml.ComponentProcessing(n)
            self._append_metabolite(cost.reactants, 'H2O', 1)
            self._append_metabolite(
                cost.products, self.default.metabolites.nmp_key(n), 1
                )
            self._append_metabolite(cost.products, 'H', 1)
            map_.component_processings.append(cost)
        return map_

    def replication_map(self):
        """
        Build replication map.

        Returns
        -------
        rba.xml.ProcessingMap
            Replication map.

        """
        map_ = rba.xml.ProcessingMap('replication')
        for n in self.default.metabolites.d_nucleotides:
            cost = rba.xml.ComponentProcessing(n)
            self._append_metabolite(
                cost.reactants, self.default.metabolites.dntp_key(n), 1
                )
            self._append_metabolite(cost.reactants, 'H2O', 1)
            self._append_metabolite(cost.products, 'PPi', 1)
            self._append_metabolite(cost.products, 'H', 1)
            map_.component_processings.append(cost)
        return map_

    def _append_metabolite(self, sr_list, key, sto):
        """Append species reference only if it was mapped to an SBML id."""
        sbml_id = self._metabolites[key].sbml_id
        if sbml_id:
            sr_list.append(rba.xml.SpeciesReference(sbml_id, sto))
