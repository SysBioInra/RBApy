"""Module importing user data needed to build the model."""

# python 2/3 compatibility
from __future__ import division, print_function, absolute_import

# global imports
import os.path
import sys
import pandas
import shutil

if sys.version_info >= (3, 9):
    from importlib import resources
else:
    import importlib_resources as resources

# local imports
from rba.prerba.pipeline_parameters import (PipelineParameters , DefaultInformation)
from rba.prerba.default_data import DefaultData
from rba.prerba import sbml_data
from rba.prerba.protein_data import ProteinData
from rba.prerba.uniprot_importer import create_uniprot_if_absent
from rba.prerba.manual_annotation import (
    CuratedMetabolites, CuratedMacrocomponents
    )
from rba.prerba.fasta_parser import RbaFastaParser
from rba.prerba import protein_export
from rba.prerba.macromolecule import Macromolecule

class UserData(object):
    """Data contained in files provided by the user."""

    def __init__(self, parameter_file, verbose=False):
        """Read data stored in filed described in parameters."""
        self.verbose = verbose
        self._parameters = PipelineParameters(parameter_file).parameters
        self.set_up_input_data_directory_structure()
        self._generate_organism_type_specific_files()
        self.default = DefaultData(species_category=self._parameters["SPECIES_CATEGORY"],default_input_directory=self.default_input_directory)
        self.uniprot_downloaded=False
        self._import_sbml_data()
        self._import_uniprot_data()
        self._import_manual_annotation()
        self.default.build_default_processing_maps(metabolite_map=self.metabolite_map,metabolite_location_map=self.metabolite_location_map,default_metabolite_location=self.protein_data.map_rba_to_sbml_compartment(self.protein_data._default_location))
        self.update_location_map_with_metabolic_compartments()
        
    def set_up_input_data_directory_structure(self):
        """
        Generates subdirectories for cellular components, machinery_composition and helper_files (if not already existing) 
        and copies their READMEs from default data.
        """
        if not os.path.isfile(self.input_path("README.rst")):
            with resources.as_file(resources.files("rba.default_data.READMES").joinpath("README.rst")) as src_path:
                shutil.copy(src_path, self.input_path("README.rst"))
        if not os.path.exists(self.input_path('cellular_components')):
            print("")
            print('Directory {} not found and therefore created.'.format(self.input_path('cellular_components')))
            os.makedirs(self.input_path('cellular_components'))
            with resources.as_file(resources.files("rba.default_data.READMES").joinpath("README_cellular_components.rst")) as src_path:
                shutil.copy(src_path, self.input_path("cellular_components/README.rst"))
        if not os.path.exists(self.input_path('machinery_composition')):
            print("")
            print('Directory {} not found and therefore created.'.format(self.input_path('machinery_composition')))
            os.makedirs(self.input_path('machinery_composition'))
            with resources.as_file(resources.files("rba.default_data.READMES").joinpath("README_machinery_composition.rst")) as src_path:
                shutil.copy(src_path, self.input_path("machinery_composition/README.rst"))
        if not os.path.exists(self.input_path('helper_files')):
            print("")
            print('Directory {} not found and therefore created.'.format(self.input_path('helper_files')))
            os.makedirs(self.input_path('helper_files'))
            with resources.as_file(resources.files("rba.default_data.READMES").joinpath("README_helper_files.rst")) as src_path:
                shutil.copy(src_path, self.input_path("helper_files/README.rst"))
        print("")

    def _generate_organism_type_specific_files(self):
        """
        Copy initial information from default data.
        """
        # define directory with default information, based on species category
        if self._parameters["SPECIES_CATEGORY"] in ["EUKARYOTE","PLANT","PROKARYOTE"]:
            self.default_input_directory="default_data/{}".format(self._parameters["SPECIES_CATEGORY"])
        else:
            print("")
            print('SPECIES_CATEGORY attribute in parameter file should be one of: \n "EUKARYOTE","PROKARYOTE" or "PLANT" ("EUKARYOTE" is assumed if not specified)')
            self.default_input_directory="default_data/EUKARYOTE"

        # read default information file
        path_default_info = resources.files('rba.' + self.default_input_directory.replace("/", ".")).joinpath("default_information.in")
        self._default_information=DefaultInformation(path_default_info).parameters       
        self._average_gene_id = self._default_information['AVERAGE_GENE_ID']
        self._macromolecule_location_separator = self._default_information['MACROMOLECULE_LOCATION_SEPARATOR']
        
        # copy process defineition file from default info directory into input data directory and read it
        if not os.path.isfile(self.input_path("process_definition.tsv")):
            with resources.as_file(resources.files('rba.' + self.default_input_directory.replace("/", ".")).joinpath("default_process_definition.tsv")) as src_path:
                shutil.copy(src_path,self.input_path("process_definition.tsv"))
        process_definition_file=pandas.read_csv(self.input_path("process_definition.tsv"),sep="\t")
        
        # store process definition as list of dictionaries, representing rows in file
        self.process_definition=[process_definition_file.loc[i,:].to_dict() for i in process_definition_file.index]
        
        # copy fasta file with tRNAs from default info directory into input data directory
        if not os.path.isfile(self.input_path("cellular_components/trnas.fasta")):
            with resources.as_file(resources.files('rba.' + self.default_input_directory.replace("/", ".")).joinpath("trnas.fasta")) as src_path:                
                shutil.copy(src_path,self.input_path("cellular_components/trnas.fasta"))

        # copy fasta files with process machinery definitions from default info directory into input data directory
        for i in self.process_definition:
            if not os.path.isfile(self.input_path('machinery_composition/{}.fasta'.format(i["PROCESS"]))):
                resource_path = resources.files('rba').joinpath(self.default_input_directory, f"{i['PROCESS']}.fasta")
                with resources.as_file(resource_path) as src_path:
                    if src_path.is_file():
                        shutil.copy(src_path, self.input_path(f'machinery_composition/{i["PROCESS"]}.fasta'))
                   
        # generate DNAs and mRNAs for each genome location, defined in default data
        genome_locations=self._default_information['GENOME_LOCATIONS'].split(",")
        self.default_dnas={loc:_generate_default_macromolecule_without_sequence(id="dna{}{}".format(self._macromolecule_location_separator,loc),location=loc,origin=loc) for loc in genome_locations}
        self.default_mrnas={loc:_generate_default_macromolecule_without_sequence(id="mrna{}{}".format(self._macromolecule_location_separator,loc),location=loc,origin=loc) for loc in genome_locations}
        self.default_genome_location=self._default_information['DEFAULT_GENOME_LOCATION']

        # import default sbml rba compartment map
        resource_path = resources.files('rba').joinpath(self.default_input_directory, "default_sbml_rba_compartment_map.tsv")
        with resources.as_file(resource_path) as src_path:
            default_sbml_rba_compartment_map = pandas.read_csv(src_path, sep="\t")
        self.default_sbml_rba_compartment_map={frozenset(default_sbml_rba_compartment_map.loc[i,"SBML_COMPARTMENT"].split(",")):default_sbml_rba_compartment_map.loc[i,"RBA_COMPARTMENT"] for i in default_sbml_rba_compartment_map.index}

    def update_location_map_with_metabolic_compartments(self):
        for met_comp in self.sbml_data.metabolic_compartments:
            self.protein_data.query_location_map(key_location=met_comp,value_location=self.sbml_data._inferred_location(set(met_comp)),comment='SBML')
            #if met_comp not in self.protein_data._location_map.data:
            #    self.protein_data._location_map.append(met_comp, self.sbml_data._inferred_location(set(met_comp)), 'SBML')
        self.protein_data._location_map.update_file(sort_by=["DATA TYPE","DATA LOCATION"])
       
    def generate_mapping_rba_to_sbml_compartments(self,rba_compartment):
        map={self.protein_data._location_map.data[i][0]:i for i in self.sbml_data.metabolic_compartments if self.protein_data._location_map.data.get(i,(None,None))[1] in ['SBML']}
        return(map.get(rba_compartment,None))

    def _import_sbml_data(self):
        if self.verbose:
            print('  Importing SBML data ...', end='')
            sys.stdout.flush()
        self.sbml_data = sbml_data.SbmlData(
            self.input_path(self._parameters['SBML_FILE']),
            external_ids=self._external_ids(),
            interface_id=self._interface_ids(),
            default_sbml_rba_compartment_map=self.default_sbml_rba_compartment_map
        )
        if self.verbose:
            print(' done')

    def input_path(self, filename):
        return os.path.join(self._input_dir(), filename)

    def _input_dir(self):
        return self._parameters['INPUT_DIR']

    def _external_ids(self):
        line = self._parameters.get('EXTERNAL_COMPARTMENTS', None)
        if line is None:
            return []
        return [e.strip() for e in line.split(',')]

    def _interface_ids(self):
        line = self._parameters.get('INTERFACE_COMPARTMENTS', None)
        if line is None:
            return []
        return set(line.split(','))

    def _import_uniprot_data(self):
        if self.verbose:
            print('  Importing UniProt data ...', end='')
            sys.stdout.flush()
        self.uniprot_downloaded=create_uniprot_if_absent(self.input_path('proteome_annotation.tsv'),self._organism_id())
        self.protein_data = ProteinData(self._input_dir(),average_gene_id=self._average_gene_id,macromolecule_location_separator=self._macromolecule_location_separator)
        # here we build "enzymatic_localised_proteins" that include the localised enzymatic proteins and the location
        # deduced from the reaction
        self._retrieve_enzymatic_localised_proteins()
        self.protein_data.update_helper_files()
        if self.verbose:
            print(' done')


    def _organism_id(self):
        return self._parameters['ORGANISM_ID']


    def _retrieve_enzymatic_localised_proteins(self):
        self.enzymatic_localised_proteins = []
        for g in self._sbml_localised_enzymatic_genes():
            gene_id = g.split(self._macromolecule_location_separator)[0]
            c_id = g.split(self._macromolecule_location_separator)[1]
            protein = self.protein_data.create_protein_from_gene_id(gene_id=gene_id,protein_id=g,location=c_id,default_genome_location=self.default_genome_location,comment_for_location_map="DoNotAdd")
            if protein:
                # update protein id with localisation
                protein.id = gene_id + self._macromolecule_location_separator + c_id
                # modify protein location of uniprot with location inferred from the reaction
                protein.location = c_id
                self.enzymatic_localised_proteins.append(protein)   

    
    def _sbml_localised_enzymatic_genes(self):
        result = []
        for enzyme in self.sbml_data.enzymes:
            result+=self._infer_enzyme_composition_localised(enzyme)
        return list(set(result))

    def _infer_enzyme_composition_localised(self,enzyme):
        """
        Generate list of enzyme subunits as localised proteins.
        Location of proteins is inferred from the location of 
        the associated metabolic species. Inference is based 
        on default mapping and once it exists on locqtion_map helper file.
        """
        result=[]
        for g in enzyme.gene_association:
            if g == '':
                continue

            enzyme_location=self.protein_data.query_location_map(key_location=' , '.join(enzyme.compartments_of_metabolites),value_location=self.sbml_data._inferred_location(set(enzyme.compartments_of_metabolites)),comment='SBML')
            
            #if ' , '.join(enzyme.compartments_of_metabolites) in self.protein_data._location_map.data:
            #    enzyme_location=self.protein_data._location_map.data[' , '.join(enzyme.compartments_of_metabolites)][0]
            #else: 
            #    enzyme_location=self.sbml_data._inferred_location(set(enzyme.compartments_of_metabolites))
            #    self.protein_data._location_map.append(' , '.join(enzyme.compartments_of_metabolites), enzyme_location, 'SBML')

            result.append(str(g +self._macromolecule_location_separator+ enzyme_location))
        return result
    

    def _import_manual_annotation(self):
        if self.verbose:
            print('  Importing manual annotation ...', end='')
            sys.stdout.flush()
        known_species = self._sbml_species_ids()
        self.macrocomponents = CuratedMacrocomponents(
            self._input_dir(), known_species
        ).data
        self.metabolite_map = self._build_metabolite_map()
        self.metabolite_location_map = self.build_metabolite_location_map()
        self.trnas = self._read_trnas(self.input_path('cellular_components/trnas.fasta'))
        if self.verbose:
            print('done')

    def _sbml_species_ids(self):
        return set([s.id for s in self.sbml_data.species])

    def build_metabolite_location_map(self):
        """Mapping of metabolic species with a list of their possible locations in metabolism"""
        metabolite_isoform_map={}
        known_species = self._sbml_species_ids()
        for species in known_species:
            proto_species=species.rsplit("_",1)[0]
            if proto_species in metabolite_isoform_map:
                metabolite_isoform_map[proto_species].append(species)
            else:
                metabolite_isoform_map[proto_species]=[species]
        return(metabolite_isoform_map)

    def _build_metabolite_map(self):
        """Map internal keys for metabolites with user-defined SBML ids."""
        known_species = self._sbml_species_ids()
        curated_data = CuratedMetabolites(self._input_dir(), known_species)
        sbml_lookup = {s.split('_', 1)[1].rsplit("_",1)[0].lower(): s.rsplit("_",1)[0] for s in known_species}
        if curated_data.file_already_existed:
            use_default_metabolites=False
        else:
            use_default_metabolites=True
        for id_, name in zip(*self._internal_species_ids_and_names(from_default_metabolites=use_default_metabolites)):
            if id_ not in curated_data.data:
                # id_ not mapped in curation file: add new entry
                if name in list(self.default.metabolites.aas_3L.values()):
                    sbml_id = sbml_lookup.get("{}__L".format(name).lower(), None)
                    if sbml_id is None:
                        sbml_id = sbml_lookup.get(name.lower(), None)
                else:
                    sbml_id = sbml_lookup.get(id_.lower(), None)

                conc = self.default.metabolites.concentration.get(id_, 0)
                # DO NOT ENTER NEW THINGS IF FILE ALREADY EXISTED BEFORE
                #if not curated_data.file_already_existed:
                curated_data.append(id_, name, sbml_id, 1, conc)
        curated_data.update_file()
        return curated_data.data

    def _internal_species_ids_and_names(self,from_default_metabolites=True):
        if from_default_metabolites:
            keys, names = self.default.metabolites.process_metabolites()
        else:
            keys=[]
            names=[]
        cofactor_info = {}
        for cofactor in self.cofactors():
            cofactor_info.setdefault(cofactor.chebi, cofactor.name)
        keys += list(cofactor_info)
        names += list(cofactor_info.values())
        return keys, names

    def _read_trnas(self, filename):
        """Read trnas in fasta file."""
        trna_data = RbaFastaParser(filename,location_separator=self._macromolecule_location_separator,metabolic_species=[sbml_species.id for sbml_species in self.sbml_data.species]).rnas
        # replace ids with user ids
        for rna in trna_data:
            user_metabolite = self.metabolite_map.get(rna.id.upper())
            if user_metabolite and user_metabolite.sbml_id:
                rna.id = user_metabolite.sbml_id
        # fuse all trnas that have the same id
        rnas = {}
        for rna in trna_data:
            aggregate_rna = rnas.get(rna.id)
            if aggregate_rna:
                aggregate_rna.sequence.append(rna.sequence)
            else:
                rna.sequence = [rna.sequence]
                rnas[rna.id] = rna
        return rnas.values()

    def average_protein_composition(self):
        # we remove non-standard amino acids from the average composition
        average_protein = self.protein_data.average_composition()
        return {aa: sto for aa, sto in average_protein.items()
                if aa in self.default.metabolites.aas}

    def average_protein_length(self):
        return sum(self.average_protein_composition().values())

    def transport_enzymes(self):
        return (e for e in self.sbml_data.enzymes if e.is_transporter)

    def metabolite_targets(self):
        result = [(id_, m.sbml_id, m.concentration) for id_, m in self.metabolite_map.items() if m.sbml_id and m.concentration]
        result += [(id_, id_, conc) for id_, conc in self.macrocomponents.items()]
        return result

    def output_dir(self):
        return self._parameters['OUTPUT_DIR']

    def export_proteins(self, filename):
        protein_export.export_proteins(
            self.input_path(filename), self.enzymatic_localised_proteins
        )

    def sbml_species(self):
        return self.sbml_data.species

    def sbml_reactions(self):
        return self.sbml_data.reactions

    def sbml_enzymes(self):
        """
        Link metabolic network and proteins.
        
        """
        enzymes = self.sbml_data.enzymes
        subunits={}
        enzyme_associated_protein_ids=[]
        removed_isoenzyme_rxns=[]
        #Check if any disabled protein (curated in protein_curation.tsv)
        #is involved with (iso)enzyme. If so remove isoenzyme and associated isrxn
        for enzyme in enzymes:
            built_enzyme=self._build_enzyme_composition_localised(enzyme)
            enzyme.composition = built_enzyme['Composition']
            subunits[enzyme.reaction]=built_enzyme['Subunits']
            if not built_enzyme['Disabled']:
                enzyme_associated_protein_ids+=[i[0] for i in enzyme.composition]
            else:
                if enzyme.reaction not in self.sbml_data.spontanous_reactions:
                    removed_isoenzyme_rxns.append(enzyme.reaction)
        
        removed_isoenzyme_rxns.sort()
        if len(removed_isoenzyme_rxns)!=0:
            for removed_isoenzyme_rxn in removed_isoenzyme_rxns:
                print("Isoenzyme removed: {}_enzyme".format(removed_isoenzyme_rxn))
                for enzyme in enzymes:
                    if enzyme.reaction==removed_isoenzyme_rxn:
                        enzymes.remove(enzyme)
                for enzyme in self.sbml_data.enzymes:
                    if enzyme.reaction==removed_isoenzyme_rxn:
                        self.sbml_data.enzymes.remove(enzyme)
                for rxn in self.sbml_data.reactions:
                    if rxn.id == removed_isoenzyme_rxn:
                        self.sbml_data.reactions.remove(rxn)

        #Generate proteins, which are enzyme subunits.#
        enzyme_associated_proteins=[]
        for protein_id in list(set(enzyme_associated_protein_ids)):
            gene_id = protein_id.split(self._macromolecule_location_separator)[0]
            loc_id = protein_id.split(self._macromolecule_location_separator)[1]
            protein = self.protein_data.create_protein_from_gene_id(gene_id=gene_id,protein_id=protein_id,location=loc_id,default_genome_location=self.default_genome_location,comment_for_location_map="DoNotAdd")
            if protein:
                protein.id = protein_id
                protein.location = loc_id
                enzyme_associated_proteins.append(protein)
        self.enzymatic_localised_proteins=enzyme_associated_proteins
        self.protein_data._curated_proteins.update_file()

        return({"Enzymes":enzymes,"Subunits":subunits})

    def _build_enzyme_composition_localised(self,enzyme):
        """
        Return enzyme composition based on localised proteins and read 
        and update manual curation file protein_curation.tsv
        """
        # Reformulate such that ifo is taken from and written to new helper file.
        gene_protein_map_enzyme={i.split(self._macromolecule_location_separator)[0]:i for i in self._infer_enzyme_composition_localised(enzyme)}
        result = {}
        subunits = {}
        isoenzyme_disabled=False
        gene_association_index=0
        for gene in enzyme.gene_association:
            gene_association_index+=1
            reference = self.protein_data.reference(gene) 
            #reference (Protein,Stoichiometry) not complex specific, since stoichiomerty is defined for protein, 
            #ergo same stoichiometry for all functions of protein.
            #Define stoichiomerty complex specifc in new curation file and use from there?
            if reference:
                gene_id=reference[0]
                enzymatic_reaction_id=enzyme.reaction.split('_duplicate_')[0]
                # here we modify the ids of the gene to include the protein localisation
                # reference is a tuple (id, stoichiometry) 
                # new referene is (localised_id,stoichiometry) 
                # but only for non-average proteins
                if self.protein_data._is_average_gene_id(gene_id):
                    #protein_id=str(gene_id.split(self._macromolecule_location_separator)[0])+self._macromolecule_location_separator+str(enzyme.location_from_reaction)
                    protein_id=str(gene_id.split(self._macromolecule_location_separator)[0])+self._macromolecule_location_separator+self.protein_data._location_map.data.get(' , '.join(enzyme.compartments_of_metabolites),(str(enzyme.location_from_reaction),None))[0]
                    
                    if protein_id not in result.keys():
                        result[protein_id]=reference[1]
                        subunits[protein_id]=gene_association_index
                    else:
                        result[protein_id]+=reference[1]

                else:
                    #check if gene-rxn (no isorxns/isoenzymes considered) combo already exists in helper file
                    if not ((self.protein_data._curated_proteins._raw_data.data["GENE"]==gene_id)&(self.protein_data._curated_proteins._raw_data.data["REACTION"]==enzymatic_reaction_id)).any():
                        # If not add it to helper file
                        self.protein_data._curated_proteins.append(gene=gene_id, 
                                                                   reaction=enzymatic_reaction_id, 
                                                                    enzyme=enzyme.reaction + '_enzyme',
                                                                    location=gene_protein_map_enzyme[gene_id].split(self._macromolecule_location_separator)[1],  
                                                                    subunit="SU_{}".format(gene_association_index),
                                                                    loc_annot=None,
                                                                    loc_sbml=None,
                                                                    stoichiometry=reference[1],
                                                                    stoich_annot_raw=self.protein_data._uniprot.line(self.protein_data._uniprot.entry(gene_id))['Subunit structure'],
                                                                    loc_annot_raw=self.protein_data.annotation_raw_location(gene=gene_id),
                                                                    id_annot=self.protein_data._uniprot.entry(gene_id)
                                                                    )
                    else:
                        # append list of isoenzymes for rxn-gene combo in helper file
                        original_enzyme_str=self.protein_data.protein_curation_get_entry(gene=gene_id,
                                                                                                   reaction=enzymatic_reaction_id,
                                                                                                   column='ISOENZYMES')
                        self.protein_data.protein_curation_add_entry(gene=gene_id,
                                                                             reaction=enzymatic_reaction_id,
                                                                             column='ISOENZYMES',
                                                                             entry=" , ".join(sorted(list(set(list(original_enzyme_str.split(" , ")+[enzyme.reaction + '_enzyme']))))))
                    
                    # Add inferred locations from different sources (metabolite locations or protein data such as uniprot)
                    self.protein_data.protein_curation_add_entry(gene=gene_id,
                                                                         reaction=enzymatic_reaction_id,
                                                                         column='LOCATION OF SBML REACTION',
                                                                         entry=gene_protein_map_enzyme[gene_id].split(self._macromolecule_location_separator)[1])
                    self.protein_data.protein_curation_add_entry(gene=gene_id,
                                                                         reaction=enzymatic_reaction_id,
                                                                         column='LOCATION FROM ANNOTATION',
                                                                         entry=" , ".join([i for i in self.protein_data.list_of_mapped_uniprot_compartments(gene=gene_id) if type(i) is str]))
    
                    # compare different inferred locations and add respective location review score as cue for manual curation in protein_curation.tsv helper file.
                    # Score 0: No conflicts
                    # Score 1: The locations inferred from different sources are unequal
                    # Score 2: The effective protein location is not equal to location inferred from metabolite locations and not among locations inferred from annotation
                    # Score 3: The effective protein location is not equal to location inferred from metabolite locations and not among locations inferred from annotation and the two information sources counterdict each other

                    listed_locations_from_annotation=self.protein_data.protein_curation_get_entry(gene=gene_id,
                                                                                                            reaction=enzymatic_reaction_id,
                                                                                                            column='LOCATION FROM ANNOTATION').split(" , ")
                    location_from_sbml=self.protein_data.protein_curation_get_entry(gene=gene_id,
                                                                                              reaction=enzymatic_reaction_id,
                                                                                              column='LOCATION OF SBML REACTION')

                    rba_location=self.protein_data.protein_curation_get_entry(gene=gene_id,
                                                                                        reaction=enzymatic_reaction_id,
                                                                                        column='RBA LOCATION')


                    if rba_location == location_from_sbml:
                        if rba_location in listed_locations_from_annotation:
                            location_review_score=0
                        else:
                            location_review_score=1
                    else:
                        if rba_location in listed_locations_from_annotation:
                            location_review_score=1
                        else:
                            if location_from_sbml in listed_locations_from_annotation:
                                location_review_score=2
                            else:
                                location_review_score=3

                    self.protein_data.protein_curation_add_entry(gene=gene_id,
                                                                            reaction=enzymatic_reaction_id,
                                                                            column='LOCATION REVIEW CODE',
                                                                            entry=location_review_score)

                    
                    # compare different inferred stoichiometries and add respective stoichiometry review score as cue for manual curation in protein_curation.tsv helper file.
                    # Score 0: Stoichiometry information was found in annotation and effective stoichiometry is equal to it
                    # Score 1: Stoichiometry information was found in annotation but effective stoichiometry is different (due to manual curation)
                    # Score 2: Stoichiometry information was not found in annotation and user defined value used as effective stoichiometry (due to manual curation)
                    # Score 3: Stoichiometry information was not found in annotation and default value (1) is used as effective stoichiometry
                    stoichiometry_from_annotation=self.protein_data._uniprot_subunits(uniprot_line=self.protein_data._uniprot.line(self.protein_data._uniprot.entry(gene_id)),default_value=None)
                    stoichiometry_rba=self.protein_data.protein_curation_get_entry(gene=gene_id,
                                                                                        reaction=enzymatic_reaction_id,
                                                                                        column='STOICHIOMETRY')
                    if stoichiometry_from_annotation: 
                        if stoichiometry_from_annotation==stoichiometry_rba:
                            stoichiometry_review_score=0
                        else:
                            stoichiometry_review_score=1
                    else:
                        if stoichiometry_rba!=1:
                            stoichiometry_review_score=2
                        else:
                            stoichiometry_review_score=3

                    self.protein_data.protein_curation_add_entry(gene=gene_id,
                                                                            reaction=enzymatic_reaction_id,
                                                                            column='STOICHIOMETRY REVIEW CODE',
                                                                            entry=stoichiometry_review_score)

                    self.protein_data._curated_proteins._raw_data._data_added=True
                        
                    #check if gene is disabled for rxn
                    if self.protein_data.protein_curation_get_entry(gene=gene_id,reaction=enzymatic_reaction_id,column='DISABLED')!=1:
                        # if no add protein to rxn-associated proteins
                        protein_id=str(gene_id)+self._macromolecule_location_separator+rba_location

                        if protein_id not in result.keys():
                            result[protein_id]=self.protein_data.protein_curation_get_entry(gene=gene_id,reaction=enzymatic_reaction_id,column='STOICHIOMETRY')
                            subunits[protein_id]=gene_association_index
                        else:
                            result[protein_id]+=self.protein_data.protein_curation_get_entry(gene=gene_id,reaction=enzymatic_reaction_id,column='STOICHIOMETRY')
                    else:
                        isoenzyme_disabled=True
        return({'Composition':list(result.items()),'Disabled':isoenzyme_disabled,'Subunits':subunits})
    
    def external_prefixes(self):
        return self.sbml_data.external_prefixes

    def compartments(self):
        return self.protein_data.compartments()

    def compartment(self, id_):
        return self.protein_data.compartment(id_)

    def cofactors(self):
        """List of all protein cofactors."""
        cofactors = []
        known_ids = set()
        for protein in self.enzymatic_localised_proteins:
            for c in protein.cofactors:
                if c.chebi not in known_ids:
                    cofactors.append(c)
                    known_ids.add(c.chebi)
        for entry in self.protein_data._cofactors.data.keys():
            for cofactor in self.protein_data._cofactors.data[entry]:
                if cofactor.chebi not in known_ids:
                    cofactors.append(cofactor)
                    known_ids.add(cofactor.chebi)
        return cofactors

def _generate_default_macromolecule_without_sequence(id,location,origin):
    mm=Macromolecule()
    mm.id = id
    mm.location = location
    mm.stoichiometry = 1
    mm.origin = origin
    return(mm)
    
