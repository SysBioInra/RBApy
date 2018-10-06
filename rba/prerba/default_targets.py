"""Module defining DefaultTargets class."""

# python 2/3 compatibility
from __future__ import division, print_function, absolute_import

# local imports
import rba.xml


class DefaultTargets(object):
    """
    Class initializing default target structures used by RBA.

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
            Default data used to initialize targets.
        metabolite_map: dict
            Map from internal metabolite ids to user-defined ids.

        """
        self.default = default_data
        self._metabolites = metabolite_map

    def translation(self, compartments):
        """
        Build translation targets.

        Parameters
        ----------
        compartments : list
            Compartment identifiers.

        Returns
        -------
        rba.xml.TargetGroup
            Translation target group.

        """
        targets = rba.xml.TargetGroup('translation_targets')
        for cpt in compartments:
            prot_id = self.default.metabolites.average_protein_id(cpt)
            target = rba.xml.TargetSpecies(prot_id)
            target.value = 'nonenzymatic_proteins_' + cpt
            targets.concentrations.append(target)
        return targets

    def transcription(self):
        """
        Build transcription targets.

        Returns
        -------
        rba.xml.TargetGroup
            Transcription targets.

        """
        targets = rba.xml.TargetGroup('transcription_targets')
        default_metabolites = self.default.metabolites
        # mrna
        target = rba.xml.TargetSpecies(default_metabolites.mrna)
        target.value = 'mrna_concentration'
        targets.concentrations.append(target)
        target = rba.xml.TargetSpecies(default_metabolites.mrna)
        target.value = 'mrna_degradation_flux'
        targets.production_fluxes.append(target)
        return targets

    def replication(self):
        """
        Build replication targets.

        Returns
        -------
        rba.xml.TargetGroup
            Replication targets.

        """
        targets = rba.xml.TargetGroup('replication_targets')
        target = rba.xml.TargetSpecies(self.default.metabolites.dna)
        target.value = 'dna_concentration'
        targets.concentrations.append(target)
        return targets

    def rna_degradation(self):
        """
        Build RNA degradation targets.

        Returns
        -------
        rba.xml.TargetGroup
            RNA degradation targets.

        """
        targets = rba.xml.TargetGroup('rna_degradation')
        target = rba.xml.TargetSpecies(self.default.metabolites.mrna)
        target.value = 'mrna_degradation_flux'
        targets.degradation_fluxes.append(target)
        return targets

    def metabolite_production(self):
        """
        Build metabolite production targets.

        Returns
        -------
        rba.xml.TargetGroup
            Metabolite production targets.

        """
        targets = rba.xml.TargetGroup('metabolite_production')
        for id_, metabolite in self._metabolites.items():
            if metabolite.sbml_id and metabolite.concentration:
                target = rba.xml.TargetSpecies(metabolite.sbml_id)
                target.value = (self.default.parameters
                                .metabolite_concentration(id_))
                targets.concentrations.append(target)
        return targets

    def macrocomponents(self, macro_fluxes):
        """
        Build macrocomponent production targets.

        Parameters
        ----------
        macro_flux : dict
            Map from molecule name to production flux.

        Returns
        -------
        rba.xml.TargetGroup
            Macrocomponent production targets.

        """
        targets = rba.xml.TargetGroup('macrocomponent_production')
        for id_ in macro_fluxes:
            target = rba.xml.TargetSpecies(id_)
            target.value = (self.default.parameters
                            .metabolite_concentration(id_))
            targets.concentrations.append(target)
        return targets

    def maintenance_atp(self, reaction_name):
        """
        Build maintenance ATP target.

        Parameters
        ----------
        reaction_name : str
            Name of maintenance ATP reaction.

        Returns
        -------
        rba.xml.TargetGroup
            Maintenance ATP targets.

        """
        targets = rba.xml.TargetGroup('maintenance_atp_target')
        target = rba.xml.TargetReaction(reaction_name)
        target.lower_bound = 'maintenance_atp'
        targets.reaction_fluxes.append(target)
        return targets
