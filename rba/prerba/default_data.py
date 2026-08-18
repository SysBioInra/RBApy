"""Module storing data and data-related functions."""

# python 2/3 compatibility
from __future__ import division, print_function, absolute_import

# local imports
import rba.xml
import pandas
import sys

if sys.version_info >= (3, 9):
    from importlib import resources
else:
    import importlib_resources as resources

from importlib_resources import files, as_file
from rba.xml import Function

GROWTH_RATE = 'growth_rate'

def build_aggregate(id_, fn_refs,aggregate_type='multiplication'):
    """Build aggregate with given identifiers and function references."""
    result = rba.xml.Aggregate(id_, aggregate_type)
    for ref in fn_refs:
        result.function_references.append(rba.xml.FunctionReference(ref))
    return result


class DefaultData(object):
    """Class holding default RBA data."""

    def __init__(self,species_category="PROKARYOTE",default_input_directory=None):
        self.default_input_directory=default_input_directory
        self.parameters = DefaultParameters(default_input_directory=self.default_input_directory)
        self.metabolites = DefaultMetabolites(default_input_directory=self.default_input_directory,species_category=species_category)
        self.atpm_reaction = 'R_maintenance_atp'


    def build_default_processing_maps(self,metabolite_map,metabolite_location_map,default_metabolite_location):
        self.processing_maps = DefaultProcessingMaps(default_metabolites=self.metabolites,
                                                     metabolite_map=metabolite_map,
                                                     metabolite_location_map=metabolite_location_map,
                                                     default_metabolite_location=default_metabolite_location)


class DefaultParameters(object):
    """Class holding default RBA parameter data."""

    def __init__(self,default_input_directory=None):
        self.default_medium_concentration=10
        self.default_transporter_efficiency='default_transporter_efficiency'
        self.default_enzyme_efficiency='default_enzyme_efficiency'

        #  Import default parameter data
        if default_input_directory is not None:
            resource_path = resources.files('rba').joinpath(default_input_directory, "default_parameters.xml")
            with resources.as_file(resource_path) as src_path:
                with open(src_path, 'r') as f:
                    self.default_parameters = rba.xml.RbaParameters().from_file(f)
                    
    def protein_fraction_functions(self, compartments):
        """
        Build functions related to relative compartment sizes. Usually values are taken from default data 
        if model compartment can be found among it. For all model compartments not found in default data, 
        a default value is calculated (assuming an equal distribution for compartments not found in default data)
        """
        default_function_ids=[i.id for i in self.default_parameters.functions._elements]
        default_aggregate_ids=[i.id for i in self.default_parameters.aggregates._elements]

        protein_fraction_ids_assigned=[]
        protein_fraction_ids_unassigned=[]
        total_fraction_assigned=0.0
        for comp in compartments:
            protein_fraction_id=self.protein_fraction_id(compartment_id=comp)
            if protein_fraction_id in default_function_ids:
                fraction_function=self.default_parameters.functions._elements_by_id[protein_fraction_id]
                if fraction_function.type == 'constant':
                    total_fraction_assigned+=float(fraction_function.parameters._elements_by_id['CONSTANT'].value)
                protein_fraction_ids_assigned.append(protein_fraction_id)
            elif protein_fraction_id in default_aggregate_ids:
                protein_fraction_ids_assigned.append(protein_fraction_id)
            else:
                protein_fraction_ids_unassigned.append(protein_fraction_id)
        
        default_protein_fraction=(1.0-total_fraction_assigned)/len(protein_fraction_ids_unassigned)
        protein_fraction_functions=[]
        for comp in compartments:
            protein_fraction_id=self.protein_fraction_id(compartment_id=comp)
            if protein_fraction_id not in list(default_function_ids+default_aggregate_ids):
                protein_fraction_functions.append(Function(protein_fraction_id, 'constant',{'CONSTANT': float(default_protein_fraction)}))
        return(protein_fraction_functions)

    def protein_pg_fraction_functions(self, compartments):
        """
        Build functions related to pg fractions sizes. Usually values are taken from default data 
        if model compartment can be found among it. For all model compartments not found explicitely in default data, 
        a default value is used)
        """
        default_function_ids=[i.id for i in self.default_parameters.functions._elements]
        default_aggregate_ids=[i.id for i in self.default_parameters.aggregates._elements]
        protein_pg_fraction_functions=[]
        for comp in compartments:
            protein_pg_fraction_id=self.non_enzymatic_fraction_id(compartment_id=comp)
            if protein_pg_fraction_id not in list(default_function_ids+default_aggregate_ids):
                if "fraction_non_enzymatic_protein__DEFAULT" in list(default_function_ids):
                    default_function=self.default_parameters.functions._elements_by_id["fraction_non_enzymatic_protein__DEFAULT"]
                    protein_pg_fraction_functions.append(Function(protein_pg_fraction_id, default_function.type,{p.id : p.value for p in default_function.parameters}))
        return(protein_pg_fraction_functions)

    def density_aggregates(self, compartments):
        """
        Build aggregate parameters related to compartment densities, serving as limits for compartment sizes. 
        """
        default_function_ids=[i.id for i in self.default_parameters.functions._elements]
        default_aggregate_ids=[i.id for i in self.default_parameters.aggregates._elements]
        density_aggregates=[]
        for comp in compartments:
            density_id="{}_density".format(comp)
            if density_id not in list(default_function_ids+default_aggregate_ids):
                density_aggregates.append(build_aggregate(density_id,
                    ['amino_acid_concentration', self.protein_fraction_id(compartment_id=comp)],
                    aggregate_type='multiplication'))
        return(density_aggregates)

    def pg_protein_aggregates(self, compartments):
        """
        Build aggregate parameters related to absolute pg amounts, serving as values for concentration targets. 
        """
        default_function_ids=[i.id for i in self.default_parameters.functions._elements]
        default_aggregate_ids=[i.id for i in self.default_parameters.aggregates._elements]
        pg_aggregates=[]
        for comp in compartments:
            pg_id='nonenzymatic_proteins_{}'.format(comp)
            if pg_id not in list(default_function_ids+default_aggregate_ids):
                pg_aggregates.append(build_aggregate('nonenzymatic_proteins_{}'.format(comp),
                    ['amino_acid_concentration',
                    'inverse_average_protein_length',
                    self.protein_fraction_id(compartment_id=comp),
                    self.non_enzymatic_fraction_id(compartment_id=comp)],
                    aggregate_type='multiplication'))
        return(pg_aggregates)

    def metabolite_concentration(self, id_):
        return id_ + '_concentration'

    def metabolite_concentration_function(self, id_, concentration):
        return Function(self.metabolite_concentration(id_),'constant', {'CONSTANT': concentration}) #ok

    def constant_function_from_id_and_value(self, id_, constant,variable=None):
        return Function(id_,'constant', {'CONSTANT': constant},variable=variable) #ok

    def default_functions(self):
        return([i for i in self.default_parameters.functions._elements if not i.id.endswith("__DEFAULT")])

    def default_aggregates(self):
        return([i for i in self.default_parameters.aggregates._elements if not i.id.endswith("__DEFAULT")])

    def transport_aggregate_id(self, reaction):
        return '{}_efficiency'.format(reaction)

    def transport_functions(self, reaction, metabolites):
        default_transport_factor_function=self.default_parameters.functions.get_by_id("transport_factor__DEFAULT")
        default_km=default_transport_factor_function.parameters.get_by_id('Km').value
        default_vmax=default_transport_factor_function.parameters.get_by_id('kmax').value
        return [
            Function(self._transport_fn_id(reaction, met),
                     'michaelisMenten', {'Km': default_km, 'kmax': default_vmax}, met)
            for met in metabolites
        ]

    def _transport_fn_id(self, reaction, metabolite):
        return '{}_{}_transport_factor'.format(reaction, metabolite)

    def transport_aggregate(self, reaction, metabolites):
        agg_fns = ([self.default_transporter_efficiency]
                   + [self._transport_fn_id(reaction, m) for m in metabolites])
        return build_aggregate(self.transport_aggregate_id(reaction), agg_fns)

    @staticmethod
    def inverse_average_protein_length(length):
        """
        Return xml structure containing average protein length.
        """
        return Function('inverse_average_protein_length', 'constant',{'CONSTANT': 1.0/length})

    @staticmethod
    def protein_fraction_id(compartment_id):
        """
        Return function identifier for protein fraction in given compartment.
        """
        return 'fraction_protein_' + compartment_id

    @staticmethod
    def non_enzymatic_fraction_id(compartment_id):
        """
        Return function identifier for non-enzymatic fraction.
        """
        return 'fraction_non_enzymatic_protein_' + compartment_id


class DefaultMetabolites(object):
    """
    Class holding default metabolite data.
    """

    def __init__(self,default_input_directory,species_category="PROKARYOTE"):

        self.species_category=species_category
        
        # import default metabolites if exist
        if default_input_directory is not None:
            resource_path = resources.files('rba').joinpath(default_input_directory, "default_metabolites.tsv")
            with resources.as_file(resource_path) as src_path:
                default_metabolites_from_file = pandas.read_csv(src_path, sep="\t")
            self.key_metabolites=dict(zip(list(default_metabolites_from_file["ID"]),list(default_metabolites_from_file["NAME"])))
            self.concentration=dict(zip(list(default_metabolites_from_file["ID"]),list(default_metabolites_from_file["CONCENTRATION"])))
        else:
            self.key_metabolites={}
            self.concentration={}

        # relative compositions of mrna and dna
        try:
            # import default dna composition 
            resource_path = resources.files('rba').joinpath(default_input_directory, "default_dna.xml")
            with resources.as_file(resource_path) as src_path:
                with open(src_path, 'r') as f:
                    self.relative_composition_dna = {i.component: i.stoichiometry for i in rba.xml.RbaDNA().from_file(f).macromolecules._elements_by_id["dna__DEFAULT"].composition._elements}
        except:
            self.relative_composition_dna = {'A': 0.2818, 'C': 0.2181, 'G': 0.2171, 'T': 0.283}
        try:
            # import default mrna composition 
            resource_path = resources.files('rba').joinpath(default_input_directory, "default_mrnas.xml")
            with resources.as_file(resource_path) as src_path:
                with open(src_path, 'r') as f:
                    self.relative_composition_mrna={i.component: i.stoichiometry for i in rba.xml.RbaRNAs().from_file(f).macromolecules._elements_by_id["mrna__DEFAULT"].composition._elements}
        except:            
            self.relative_composition_mrna = {'A': 0.2818, 'C': 0.2181, 'G': 0.2171, 'U': 0.283}

        # key metabolites and components
        self.aas_3L = {
            'A': 'ala', 'C': 'cys', 'D': 'asp', 'E': 'glu', 'F': 'phe',
            'G': 'gly', 'H': 'his', 'I': 'ile', 'K': 'lys', 'L': 'leu',
            'fM': 'fmet', 'M': 'met', 'N': 'asn', 'P': 'pro', 'Q': 'gln',
            'R': 'arg', 'S': 'ser', 'T': 'thr', 'V': 'val', 'W': 'trp',
            'Y': 'tyr'
            }
        
        self.aas = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
        self.aa_fM = 'fM'
        self.nucleotides = ['A', 'C', 'G', 'U']
        self.d_nucleotides = ['A', 'C', 'G', 'T']

        ###############

    def charged_trna_key(self, aa):
        """
        Return internal identifier of charged trna.
        """
        return self.aas_3L[aa].upper() + 'TRNA'

    def charged_trna_name(self, aa):
        """
        Return name of charged trna.
        """
        return 'Charged trna ' + self.aas_3L[aa]

    def uncharged_trna_key(self, aa):
        """
        Return internal identifier of uncharged trna.
        """
        return 'TRNA' + self.aas_3L[aa].upper()

    def uncharged_trna_name(self, aa):
        """
        Return name of uncharged trna.
        """
        return 'Uncharged trna ' + self.aas_3L[aa]

    @staticmethod
    def ntp_key(letter):
        """
        Return internal identifier of triphosphote nucleotide.
        """
        return letter + 'TP'

    @staticmethod
    def ndp_key(letter):
        """
        Return internal identifier of diphosphate nucleotide.
        """
        return letter + 'DP'

    @staticmethod
    def nmp_key(letter):
        """
        Return internal identifier of monophosphate nucleotide.
        """
        return letter + 'MP'
        
    @staticmethod
    def dntp_key(letter):
        """
        Return internal identifier of triphosphate deoxynucleotide.
        """
        return 'd' + letter + 'TP'

    def process_metabolites(self):
        """
        Return internal ids and names of metabolites involved in processes.
        """
        # metabolites to retrieve
        keys = []  # internal id of metabolite
        names = []  # name of metabolite
        # methionine
        if self.species_category=="PROKARYOTE":
            keys.append('MET')
            names.append('Methionine')

        # free amino acids
        keys += ["{} ({})".format(aa,self.aas_3L[aa]) for aa in self.aas]
        names += [self.aas_3L[aa] for aa in self.aas]
        # charged + uncharged trnas
        names += [self.charged_trna_name(aa) for aa in self.aas]
        keys += [self.charged_trna_key(aa) for aa in self.aas]
        names += [self.uncharged_trna_name(aa) for aa in self.aas]
        keys += [self.uncharged_trna_key(aa) for aa in self.aas]
        keys.append(self.charged_trna_key(self.aa_fM))
        names.append(self.charged_trna_name(self.aa_fM))
        # nucleotides
        keys += [self.ntp_key(n) for n in self.nucleotides]
        names += [self.ntp_key(n) for n in self.nucleotides]
        keys += [self.ndp_key(n) for n in self.nucleotides]
        names += [self.ndp_key(n) for n in self.nucleotides]
        keys += [self.nmp_key(n) for n in self.nucleotides]
        names += [self.nmp_key(n) for n in self.nucleotides]
        keys += [self.dntp_key(n) for n in self.d_nucleotides]
        names += [self.dntp_key(n) for n in self.d_nucleotides]
        # key metabolites
        for met_id, name in self.key_metabolites.items():
            keys.append(met_id)
            names.append(name)
        return keys, names


class DefaultProcessingMaps(object):
    """Class holding default processing maps, for key cellular processes."""

    def __init__(self,default_metabolites,metabolite_map,metabolite_location_map,default_metabolite_location):
        self.default_metabolites=default_metabolites
        self._metabolites = metabolite_map
        self._metabolite_location_map = metabolite_location_map
        self.default_metabolite_location=default_metabolite_location

    def default_map(self,map_id,macromolecule_type,non_spontaneous=True):
        """
        Build default map with 1 machinery cost per component.

        Returns
        -------
        rba.xml.ProcessingMap
            Folding map.

        """
        if macromolecule_type =="protein":
            components=self.default_metabolites.aas
        elif macromolecule_type =="rna":
            components=self.default_metabolites.nucleotides
        elif macromolecule_type =="dna":
            components=self.default_metabolites.d_nucleotides
        map_ = rba.xml.ProcessingMap(map_id)
        if non_spontaneous:
            machinery_cost=1
        else:
            machinery_cost=0
        for component in components:
            component_processing=rba.xml.ComponentProcessing(component, machinery_cost)
            map_.component_processings.append(component_processing)
        return map_


    def default_map_atp_dependent(self,map_id,processing_location,macromolecule_type,non_spontaneous=True):
        """
        Build default map with 1 machinery cost and one atp turnover per component.

        Returns
        -------
        rba.xml.ProcessingMap
            Folding map.

        """
        if macromolecule_type =="protein":
            components=self.default_metabolites.aas
        elif macromolecule_type =="rna":
            components=self.default_metabolites.nucleotides
        elif macromolecule_type =="dna":
            components=self.default_metabolites.d_nucleotides
        map_ = rba.xml.ProcessingMap(map_id)
        if non_spontaneous:
            machinery_cost=1
        else:
            machinery_cost=0
        for component in components:
            component_processing=rba.xml.ComponentProcessing(component, machinery_cost)
            self._append_metabolite(component_processing.reactants, 'ATP', 1.0, target_location=processing_location)
            self._append_metabolite(component_processing.reactants, 'H2O', 1.0, target_location=processing_location)
            self._append_metabolite(component_processing.products, 'ADP', 1.0, target_location=processing_location)
            self._append_metabolite(component_processing.products, 'Pi', 1.0, target_location=processing_location)
            map_.component_processings.append(component_processing)
        return map_


    def default_translation_map(self,map_id,processing_location, cofactors, is_prokaryote=True,non_spontaneous=True):
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
        map_ = rba.xml.ProcessingMap(map_id)

        if non_spontaneous:
            machinery_cost=1
        else:
            machinery_cost=0
        if is_prokaryote:
            # constant costs
            reactants = map_.constant_processing.reactants
            self._append_metabolite(reactants, self.default_metabolites.charged_trna_key('fM'), 1, target_location=processing_location)
            self._append_metabolite(reactants, 'GTP', 1, target_location=processing_location)
            self._append_metabolite(reactants, 'H2O', 2, target_location=processing_location)
            products = map_.constant_processing.products
            self._append_metabolite(products,self.default_metabolites.uncharged_trna_key('M'), 1, target_location=processing_location)
            self._append_metabolite(products, 'MET', 1, target_location=processing_location)
            self._append_metabolite(products, 'FOR', 1, target_location=processing_location)
            self._append_metabolite(products, 'GDP', 1, target_location=processing_location)
            self._append_metabolite(products, 'Pi', 1, target_location=processing_location)
            self._append_metabolite(products, 'H', 1, target_location=processing_location)
        # amino acids
        for aa in self.default_metabolites.aas:
            cost = rba.xml.ComponentProcessing(aa, machinery_cost)
            self._append_metabolite(cost.reactants, self.default_metabolites.charged_trna_key(aa), 1, target_location=processing_location)
            self._append_metabolite(cost.reactants, 'GTP', 2, target_location=processing_location)
            self._append_metabolite(cost.reactants, 'H2O', 2, target_location=processing_location)
            self._append_metabolite(cost.products, self.default_metabolites.uncharged_trna_key(aa), 1, target_location=processing_location)
            self._append_metabolite(cost.products, 'GDP', 2, target_location=processing_location)
            self._append_metabolite(cost.products, 'Pi', 2, target_location=processing_location)
            self._append_metabolite(cost.products, 'H', 2, target_location=processing_location)
            map_.component_processings.append(cost)
        # cofactors
        for cofactor in cofactors:
            cost = rba.xml.ComponentProcessing(cofactor.chebi)
            self._append_metabolite(cost.reactants, cofactor.chebi, 1, target_location=processing_location)
            map_.component_processings.append(cost)
        return map_


    def default_transcription_map(self,map_id,processing_location,non_spontaneous=True):
        """
        Build transcription map.

        Returns
        -------
        rba.xml.ProcessingMap
            Transcription map.

        """            
        map_ = rba.xml.ProcessingMap(map_id)
        if non_spontaneous:
            machinery_cost=1
        else:
            machinery_cost=0
        for n in self.default_metabolites.nucleotides:
            cost = rba.xml.ComponentProcessing(n,machinery_cost)
            self._append_metabolite(cost.reactants, self.default_metabolites.ntp_key(n), 1, target_location=processing_location)
            self._append_metabolite(cost.reactants, 'H2O', 1, target_location=processing_location)
            self._append_metabolite(cost.products, 'PPi', 1, target_location=processing_location)
            self._append_metabolite(cost.products, 'H', 1, target_location=processing_location)
            map_.component_processings.append(cost)
        return map_


    def default_replication_map(self,map_id,processing_location,non_spontaneous=True):
        """
        Build replication map.

        Returns
        -------
        rba.xml.ProcessingMap
            Replication map.

        """
        map_ = rba.xml.ProcessingMap(map_id)
        if non_spontaneous:
            machinery_cost=1
        else:
            machinery_cost=0
        for n in self.default_metabolites.d_nucleotides:
            cost = rba.xml.ComponentProcessing(n,machinery_cost)
            self._append_metabolite(cost.reactants, self.default_metabolites.dntp_key(n), 1, target_location=processing_location)
            self._append_metabolite(cost.reactants, 'H2O', 1, target_location=processing_location)
            self._append_metabolite(cost.products, 'PPi', 1, target_location=processing_location)
            self._append_metabolite(cost.products, 'H', 1, target_location=processing_location)
            map_.component_processings.append(cost)
        return map_


    def default_rna_degradation_map(self,map_id,processing_location,non_spontaneous=True):
        """
        Build RNA degradation map.

        Returns
        -------
        rba.xml.ProcessingMap
            RNA degradation map.

        """
        map_ = rba.xml.ProcessingMap(map_id)
        if non_spontaneous:
            machinery_cost=1
        else:
            machinery_cost=0
        for n in self.default_metabolites.nucleotides:
            cost = rba.xml.ComponentProcessing(n,machinery_cost)
            self._append_metabolite(cost.reactants, 'H2O', 1, target_location=processing_location)
            self._append_metabolite(cost.products, self.default_metabolites.nmp_key(n), 1, target_location=processing_location)
            self._append_metabolite(cost.products, 'H', 1, target_location=processing_location)
            map_.component_processings.append(cost)
        return map_


    def default_protein_degradation_map(self,map_id,processing_location,cofactors,non_spontaneous=True):
        map_ = rba.xml.ProcessingMap(map_id)
        if non_spontaneous:
            machinery_cost=1
        else:
            machinery_cost=0
        # amino acids
        for aa in self.default_metabolites.aas:
            component_processing = rba.xml.ComponentProcessing(aa, machinery_cost)
            self._append_metabolite(component_processing.reactants, 'ATP', 1, target_location=processing_location)
            self._append_metabolite(component_processing.reactants, 'H2O', 1, target_location=processing_location)
            self._append_metabolite(component_processing.products, "{} ({})".format(aa,self.default_metabolites.aas_3L[aa]), 1, target_location=processing_location)
            self._append_metabolite(component_processing.products, 'ADP', 1, target_location=processing_location)
            self._append_metabolite(component_processing.products, 'Pi', 1, target_location=processing_location)
            map_.component_processings.append(component_processing)
        # cofactors
        for cofactor in cofactors:
            component_processing = rba.xml.ComponentProcessing(cofactor.chebi)
            self._append_metabolite(component_processing.products, cofactor.chebi, 1, target_location=processing_location)
            map_.component_processings.append(component_processing)
        return map_


    def _append_metabolite(self, sr_list, key, sto, target_location=""):
        """Append species reference only if it was mapped to an SBML id."""
        proto_sbml_id = self._metabolites[key].sbml_id
        multiplier = self._metabolites[key].stoichiometry # STOICHIOMETRY in metabolite.tsv
        if proto_sbml_id:
            located_sbml_ids=self._metabolite_location_map[proto_sbml_id]
            if "{}_{}".format(proto_sbml_id,target_location) in located_sbml_ids:
                sr_list.append(rba.xml.SpeciesReference("{}_{}".format(proto_sbml_id,target_location), sto*multiplier))
            elif "{}_{}".format(proto_sbml_id,self.default_metabolite_location) in located_sbml_ids:
                sr_list.append(rba.xml.SpeciesReference("{}_{}".format(proto_sbml_id,self.default_metabolite_location), sto*multiplier))
            else:
                sr_list.append(rba.xml.SpeciesReference(located_sbml_ids[0], sto*multiplier))

