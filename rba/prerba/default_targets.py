"""Module defining DefaultTargets class."""

# python 2/3 compatibility
from __future__ import division, print_function, absolute_import

import sys

if sys.version_info >= (3, 9):
    from importlib import resources
else:
    import importlib_resources as resources

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

    def __init__(self, user_data):
        """
        Build object from default data and metabolite map.

        Parameters
        ----------
        user_data : rba.prerba.default_data.User)
            USer data used to initialize targets.

        """
        self.user_data = user_data
        # import default target
        if self.user_data.default_input_directory is not None:
            resource_path = resources.files('rba').joinpath(self.user_data.default_input_directory, "default_targets.xml")
            with resources.as_file(resource_path) as src_path:
                with open(src_path, 'r') as f:
                    self.default_targets =rba.xml.RbaTargets().from_file(f)


    def transcription(self):
        """
        Build transcription targets.

        Returns
        -------
        rba.xml.TargetGroup
            Transcription targets.

        """
        return(self.default_targets.target_groups._elements_by_id.get('transcription_targets',None))


    def replication(self):
        """
        Build replication targets.

        Returns
        -------
        rba.xml.TargetGroup
            Replication targets.

        """
        return(self.default_targets.target_groups._elements_by_id.get('replication_targets',None))


    def rna_degradation(self):
        """
        Build RNA degradation targets.

        Returns
        -------
        rba.xml.TargetGroup
            RNA degradation targets.

        """
        return(self.default_targets.target_groups._elements_by_id.get('rna_degradation',None))


    def translation(self):
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
        for comp in self.user_data.compartments():
            prot_id = self.user_data.protein_data.average_protein_id(comp)
            target = rba.xml.TargetSpecies(prot_id)
            target.value = 'nonenzymatic_proteins_{}'.format(comp)
            targets.concentrations.append(target)
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
        internal_metabolite_species=[i.id for i in self.user_data.sbml_data.species if not i.boundary_condition]

        for id_, metabolite in self.user_data.metabolite_map.items():
            if metabolite.sbml_id and metabolite.concentration:
                # define metabolite, located in cytoplasm as target species. 
                located_internal_sbml_ids=[i for i in self.user_data.metabolite_location_map[metabolite.sbml_id] if i in internal_metabolite_species]
                target_metabolite=None
                if "{}_{}".format(metabolite.sbml_id,self.user_data.generate_mapping_rba_to_sbml_compartments(rba_compartment=self.user_data.compartment('Cytoplasm'))) in located_internal_sbml_ids:
                    target_metabolite="{}_{}".format(metabolite.sbml_id,self.user_data.generate_mapping_rba_to_sbml_compartments(rba_compartment=self.user_data.compartment('Cytoplasm')))
                #if metabolite is not present in cytoplasm, chose the one in default location if present in this copmpartmen
                elif "{}_{}".format(metabolite.sbml_id,self.user_data.generate_mapping_rba_to_sbml_compartments(rba_compartment=self.user_data.protein_data._default_location)) in located_internal_sbml_ids:
                    target_metabolite="{}_{}".format(metabolite.sbml_id,self.user_data.generate_mapping_rba_to_sbml_compartments(rba_compartment=self.user_data.compartment(self.user_data.protein_data._default_location)))
                #otherwise chose species, regardles of location
                elif len(located_internal_sbml_ids)!=0:
                    target_metabolite=located_internal_sbml_ids[0]
                if target_metabolite is not None:
                    target = rba.xml.TargetSpecies(target_metabolite)
                    target.value = (self.user_data.default.parameters.metabolite_concentration(metabolite.sbml_id))
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
            target.value = (self.user_data.default.parameters.metabolite_concentration(id_))
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
