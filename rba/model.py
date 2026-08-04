"""Module defining RbaModel class."""

# python 2/3 compatibility
from __future__ import division, print_function, absolute_import

# global imports
import os
import shutil
import sys

if sys.version_info >= (3, 9):
    from importlib import resources
else:
    import importlib_resources as resources

# local imports
import rba
from rba.utils import efficiencies

class RbaModel(object):
    """
    Class holding RBA model.

    Attributes:
        metabolism: xml structure with metabolism data.
        parameters: xml structure with parameter data.
        compartments: xml structure with compartment data.
        proteins: xml structure with protein data.
        enzymes: xml structure with enzyme data.
        rnas: xml structure with rna data.
        dna: xml structure with dna data.
        other_macromolecules: xml structure with other_macromolecules data.
        processes: xml structure with process data.
        targets: xml structure with target data.
        custom_constraints: xml structure with custom_constraints data.
        medium: dictionary mapping metabolite prefixes with their medium
            concentration.
        output_dir: path to directory where model files should be written.

    """

    def __init__(self):
        """Build object with empty structures."""
        self.metabolism = rba.xml.RbaMetabolism()
        self.compartments = rba.xml.RbaCompartments()
        self.parameters = rba.xml.RbaParameters()
        self.proteins = rba.xml.RbaProteins()
        self.enzymes = rba.xml.RbaEnzymes()
        self.rnas = rba.xml.RbaRNAs()
        self.dna = rba.xml.RbaDNA()
        self.other_macromolecules = rba.xml.RbaMacromolecules()
        self.processes = rba.xml.RbaProcesses()
        self.targets = rba.xml.RbaTargets()
        self.custom_constraints = rba.xml.RbaCustomConstraints()
        self.medium = {}
        self.output_dir = ''
        self._constraint_matrix = None
        self.get_metadata(os.path.join(self.output_dir, 'metadata.tsv'))

    @classmethod
    def from_data(cls, params_file, verbose=False):
        """
        Make object from data directory (SBML, FASTA and helper files).

        Parameters
        ----------
        params_file : str
            Path to params.in file.

        verbose : bool, optional
            Whether to display status information
        """
        builder = rba.ModelBuilder(params_file, verbose=verbose)
        builder.export_proteins('helper_files/protein_summary.tsv')
        model = builder.build_model()
    
        return model

    @classmethod
    def from_xml(cls, input_dir):
        """
        Make object from xml files.

        Parameters
        ----------
        input_dir : str
            Path to directory containing RBA XML files.

        """
        model_file_index=ModelFileIndex()
        model_file_index.read_from_file(os.path.join(input_dir, 'model_file_index.in'))
        obj = cls()
        obj.output_dir = input_dir
        obj.get_metadata(os.path.join(input_dir, 'metadata.tsv'))
        obj.parameters = rba.xml.RbaParameters().from_file(
            open(os.path.join(input_dir, model_file_index.parameters.get('parameters')))
        )
        obj.metabolism = rba.xml.RbaMetabolism().from_file(
            open(os.path.join(input_dir, model_file_index.parameters.get('metabolism')))
        )
        obj.proteins = rba.xml.RbaProteins().from_file(
            open(os.path.join(input_dir, model_file_index.parameters.get('proteins')))
        )
        obj.rnas = rba.xml.RbaRNAs().from_file(
            open(os.path.join(input_dir, model_file_index.parameters.get('rnas')))
        )
        obj.dna = rba.xml.RbaDNA().from_file(
            open(os.path.join(input_dir, model_file_index.parameters.get('dna')))
        )
        obj.other_macromolecules = rba.xml.RbaMacromolecules().from_file(
            open(os.path.join(input_dir, model_file_index.parameters.get('other_macromolecules')))
        )
        obj.processes = rba.xml.RbaProcesses().from_file(
            open(os.path.join(input_dir, model_file_index.parameters.get('processes')))
        )
        obj.targets = rba.xml.RbaTargets().from_file(
            open(os.path.join(input_dir, model_file_index.parameters.get('targets')))
        )
        obj.enzymes = rba.xml.RbaEnzymes().from_file(
            open(os.path.join(input_dir, model_file_index.parameters.get('enzymes')))
        )
        obj.custom_constraints = rba.xml.RbaCustomConstraints().from_file(
            open(os.path.join(input_dir, model_file_index.parameters.get('custom_constraints')))
        )
        obj.compartments = rba.xml.RbaCompartments().from_file(
            open(os.path.join(input_dir, model_file_index.parameters.get('compartments')))
        )
        obj.set_medium(os.path.join(input_dir, model_file_index.parameters.get('medium')))

        if os.path.isfile(os.path.join(input_dir, 'density.xml')):
            obj.density = rba.xml.RbaDensity().from_file(
                open(os.path.join(input_dir, 'density.xml'))
            )

        return obj

    def get_metadata(self, file_name):
        """
        Create meta-data structure.
        If no file exists to import meta data, it is initiated empty.

        Args:
            file_name: path to file containing metadata.
        """
        metadata = {"Date_first_model_generation":None,
                    "Date_last_model_generation":None,
                    "RBApy_version_model_generation":None,
                    "Date_download_proteome_annotation":None,
                    "Date_of_model_update_to_v3_or_higher":None,
                    "RBApy_version_model_update":None,
                    "RBAML_version":None,
                    "Model":None,
                    "Organism/Tissue":None,
                    "Taxon_ID":None,
                    "Author(s)":None,
                    "Copyright_holder":None,
                    "License":None,
                    "Localised_Macromolecules":False,
                    "Macromolecule_location_separator":None}
        if os.path.isfile(file_name):
            with open(file_name, 'r') as input_stream:
                for line in input_stream:
                    [key, value] = line.rstrip().split('\t')
                    metadata.update({key:value})
        self.meta_data = metadata

    def set_medium(self, file_name):
        """
        Set medium concentrations according to file.

        Args:
            file_name: path to file containing medium concentrations in a tab
                separated format. File is supposed to contain a header, one
                column with metabolite prefixes and one column with
                concentrations values.
        """
        concentrations = {}
        with open(file_name, 'r') as input_stream:
            # skip header
            next(input_stream)
            for line in input_stream:
                [met, conc] = line.rstrip().split('\t')[:2] #make sure only first two columns of file are used
                concentrations[met] = float(conc)
        self.medium = concentrations
        if self._constraint_matrix:
            self._constraint_matrix.set_medium(concentrations)

    def generate_mean_composition_machinery(self):
        """
        Generates mean composition enzymes from isoenzymes.
        """
        enzyme_ids=[i.id for i in list(self.enzymes.enzymes._elements_by_id.values())]

        # generating proto-isoenzyme map 
        # {'R_A_enzyme':['R_A_enzyme','R_A_duplicate_2_enzyme'],...}:
        isoenzyme_map={}
        for e_id in enzyme_ids:
            if "_duplicate_" in e_id:
                proto_id=e_id.split("_duplicate_")[0]+"_enzyme"
            else:
                proto_id=e_id
            if proto_id not in isoenzyme_map:
                isoenzyme_map[proto_id]=[e_id]
            else:
                isoenzyme_map[proto_id].append(e_id)

        ##### Generate mean-composition enzyme for each reaction and remove duplicates
        mean_composition_enzymes = rba.xml.RbaEnzymes()
        # Count the number of occurences of each reactant and product (subunits,Metabolite, etc) 
        # among all isoenzymes (respecting stoichiometry). Value, divided by numer of isoenzymes
        # is stoichiometric factor in mean-composition enzyme:
        for proto_enzyme_id in isoenzyme_map.keys():
            #Count the number of occurences of each reactant and product and count isoenzymes:
            reactants={}
            products={}
            reactant_comments={}
            product_comments={}
            isoenzyme_count=0
            for iso_enzyme_id in isoenzyme_map[proto_enzyme_id]:
                isoenzyme=self.enzymes.enzymes._elements_by_id[iso_enzyme_id]
                for reactant in isoenzyme.machinery_composition.reactants:
                    if reactant.species in reactants:
                        reactant_comments[reactant.species].append(reactant.comment)
                        reactants[reactant.species]+=reactant.stoichiometry
                    else:
                        reactant_comments[reactant.species]=[reactant.comment]
                        reactants[reactant.species]=reactant.stoichiometry
                for product in isoenzyme.machinery_composition.products:
                    if product.species in products:
                        product_comments[reactant.species].append(product.comment)
                        products[product.species]+=product.stoichiometry
                    else:
                        product_comments[reactant.species]=[product.comment]
                        products[product.species]=product.stoichiometry
                isoenzyme_count+=1

            # Generate mean-composition enzyme object:

            isoenzymes=isoenzyme_map[proto_enzyme_id]
            isoenzymes.sort()        
            mean_compo_enzyme = rba.xml.Enzyme(proto_enzyme_id,
                                               str(proto_enzyme_id.split("_enzyme")[0]), 
                                               self.enzymes.enzymes._elements_by_id[isoenzymes[0]].forward_efficiency,
                                               self.enzymes.enzymes._elements_by_id[isoenzymes[0]].backward_efficiency)
            for reactant in reactants.keys():
                if None in reactant_comments[reactant]:
                    mean_compo_enzyme.machinery_composition.reactants.append(rba.xml.SpeciesReference(species=reactant, stoichiometry=reactants[reactant]/isoenzyme_count))
                else:
                    mean_compo_enzyme.machinery_composition.reactants.append(rba.xml.SpeciesReference(species=reactant, stoichiometry=reactants[reactant]/isoenzyme_count,comment="/".join(list(set(reactant_comments[reactant])))))
            for product in products.keys():
                if None in product_comments[product]:
                    mean_compo_enzyme.machinery_composition.products.append(rba.xml.SpeciesReference(species=product, stoichiometry=products[product]/isoenzyme_count))
                else:
                    mean_compo_enzyme.machinery_composition.products.append(rba.xml.SpeciesReference(species=product, stoichiometry=products[product]/isoenzyme_count,comment="/".join(list(set(product_comments[product])))))
            # append to other mean-composition enzymes:
            mean_composition_enzymes.enzymes.append(mean_compo_enzyme)
        
        ### Generate respective metabolism object, by only consideren non-duplicated reactions:
        mean_composition_metabolism = rba.xml.RbaMetabolism()
        mean_composition_metabolism.species=self.metabolism.species
        rxns_added=[]
        for rxn in self.metabolism.reactions:
            rxn_id=rxn.id
            if rxn_id.split("_duplicate")[0] not in rxns_added:
                rxn.id=rxn_id.split("_duplicate")[0]
                mean_composition_metabolism.reactions.append(rxn)
                rxns_added.append(rxn_id.split("_duplicate")[0])

        return(mean_composition_enzymes , mean_composition_metabolism)

    def write(self, output_dir=None, generate_mean_composition_model=False):
        """
        Write rba files in XML format, model_file_index and medium and metadata as tsv.

        Parameters
        ----------
        output_dir : str
            Path to directory where files should be written. If
            specified, output_dir attribute is overriden, otherwise value
            currently stored in output_dir attribute is used.
        generate_mean_composition_model : bool
            Whether to generate xml files for enzymes and metabolism, 
            without all duplicated reactions/enzymes (isoenzymes) 
            and mean-composition enzymes. 
            Filenames are mean_composition_model_enzymes and
            mean_composition_model_metabolism respectively.

        """
        if output_dir:
            self.output_dir = output_dir

        if not os.path.isfile(self._output("README.rst")):
            with resources.as_file(resources.files("rba.default_data.READMES").joinpath("README.rst")) as src_path:
                shutil.copy(src_path, self._output("README.rst"))
        if not os.path.exists(self._output('model')):
            os.makedirs(self._output('model'))
            with resources.as_file(resources.files("rba.default_data.READMES").joinpath("README_model.rst")) as src_path:
                shutil.copy(src_path, self._output("model/README.rst"))
        if not os.path.exists(self._output('outputs')):
            os.makedirs(self._output('outputs'))
            with resources.as_file(resources.files("rba.default_data.READMES").joinpath("README_outputs.rst")) as src_path:
                shutil.copy(src_path, self._output("outputs/README.rst"))
        if not os.path.exists(self._output('other')):
            os.makedirs(self._output('other'))
            with resources.as_file(resources.files("rba.default_data.READMES").joinpath("README_other.rst")) as src_path:
                shutil.copy(src_path, self._output("other/README.rst"))

        model_file_index=ModelFileIndex()
        model_file_index.parameters={'compartments':'model/compartments.xml', 
                                     'metabolism':'model/metabolism.xml', 
                                     'enzymes':'model/enzymes.xml',
                                     'proteins':'model/proteins.xml',
                                     'rnas':'model/rnas.xml',
                                     'dna':'model/dna.xml',
                                     'other_macromolecules':'model/other_macromolecules.xml',
                                     'processes':'model/processes.xml',
                                     'targets':'model/targets.xml',
                                     'parameters':'model/parameters.xml',
                                     'custom_constraints':'model/custom_constraints.xml',
                                     'medium':'model/medium.tsv'}
        self.metabolism.write(self._output(model_file_index.parameters['metabolism']),meta_data=self.meta_data)
        self.compartments.write(self._output(model_file_index.parameters['compartments']),meta_data=self.meta_data)
        self.proteins.write(self._output(model_file_index.parameters['proteins']),meta_data=self.meta_data)
        self.rnas.write(self._output(model_file_index.parameters['rnas']),meta_data=self.meta_data)
        self.dna.write(self._output(model_file_index.parameters['dna']),meta_data=self.meta_data)
        self.other_macromolecules.write(self._output(model_file_index.parameters['other_macromolecules']),meta_data=self.meta_data)
        self.enzymes.write(self._output(model_file_index.parameters['enzymes']),meta_data=self.meta_data)
        self.parameters.write(self._output(model_file_index.parameters['parameters']),meta_data=self.meta_data)
        self.processes.write(self._output(model_file_index.parameters['processes']),meta_data=self.meta_data)
        self.targets.write(self._output(model_file_index.parameters['targets']),meta_data=self.meta_data)
        self.custom_constraints.write(self._output(model_file_index.parameters['custom_constraints']),meta_data=self.meta_data)
        
        # initial conditions (medium concentrations)
        with open(self._output(model_file_index.parameters['medium']), 'w') as output:
            output.write('Metabolite\tConcentration\tType\tUnit\tComment\n')
            for met, conc in self.medium.items():
                output.write('{}\t{}\t{}\t{}\t{}\n'.format(met, conc,"Chemical species","mmol/L","Boundary metabolite from SBML"))
            if "Temperature" not in self.medium.keys():
                output.write('{}\t{}\t{}\t{}\t{}\n'.format("Temperature", 20,"Temperature","°Celsius","Environmental condition"))

        with open(self._output('metadata.tsv'), 'w') as output:
            for key, value in self.meta_data.items():
                if value is not None:
                    output.write('{}\t{}\n'.format(key, value))
                else:
                    output.write('{}\t{}\n'.format(key, "TBD"))

        if generate_mean_composition_model:
            mean_composition_enzymes , mean_composition_metabolism = self.generate_mean_composition_machinery()
            mean_composition_enzymes.write(self._output('model/mean_composition_model_enzymes.xml'),meta_data=self.meta_data)
            mean_composition_metabolism.write(self._output('model/mean_composition_model_metabolism.xml'),meta_data=self.meta_data)
            model_file_index.parameters.update({'# metabolism':'model/mean_composition_model_metabolism.xml'})
            model_file_index.parameters.update({'# enzymes':'model/mean_composition_model_enzymes.xml'})
            
        model_file_index.write_to_file(file_path=self._output('model_file_index.in'))

    def _output(self, file_name):
        """Return full path to file contained in output directory."""
        return os.path.join(self.output_dir, file_name)

    def solve(self, recompute_matrices=True, lp_solver=None,
              mu_min=0., mu_max=2.5,
              bissection_tol=1e-6, max_bissection_iters=None,
              verbose=False,solve_grid=False):
        """
        Solve RBA model.

        Parameters
        ----------
        recompute_matrices : bool, optional
            If the model is solved several time, recompute matrices defining
            the optimality problem (True by default). This parameter should be
            set to False when only the medium composition changes (medium
            concentrations do not appear in the matrices).

        lp_solver : str, optional
            LP solver (``cplex``, ``glpk``, ``gurobi``, ``scipy`` ``swiglpk``)

        mu_min : float, optional:
            Minimum μ to check

        mu_max : float, optional:
            Maximum μ to check

        bissection_tol : float, optional
            Tolerance for bissection

        max_bissection_iters : int, optional
            Maximum number of iterations for bissection

        verbose : bool, optional
            Whether to display status information

        Returns
        -------
        rba.utils.Results object containing optimal growth rate and fluxes.

        """
        if recompute_matrices or self._constraint_matrix is None:
            self._constraint_matrix = rba.ConstraintMatrix(self)
        solver = rba.Solver(self._constraint_matrix,
                            lp_solver=lp_solver,
                            mu_min=mu_min,
                            mu_max=mu_max,
                            bissection_tol=bissection_tol,
                            max_bissection_iters=max_bissection_iters,
                            verbose=verbose)
        if solve_grid:
            solver.solve_grid()
        else:        
            solver.solve()

        return rba.Results(self, self._constraint_matrix, solver)

    def set_enzyme_efficiencies(self, file_name):
        efficiencies.set_efficiencies(self, file_name)


class ModelFileIndex(object):
    """
    Class storing model_file_index parameters.

    Attributes
    ----------
    obligatory_tags : list of str
        tags that a model_file_index file must contain.
    parameters: dict
        Dictonary mapping parameter tags with their value. Optional tags
        may be omitted.

    """

    def __init__(self):
        """
        Build model_file_index object without parameters
        """
        self.obligatory_tags = ['compartments', 'metabolism', 'enzymes',
                                'proteins','rnas','dna','other_macromolecules'
                                ,'processes','targets','parameters',
                                'custom_constraints','medium']
        self.parameters = {}

    def read_from_file(self,file_path):
        """
        Build model_file_index object from file path.

        Parameters
        ----------
        file_path: str
            path to parameter file.

        """
        try:
            with open(file_path, 'r') as input_stream:
                # parse file
                for line in input_stream:
                    line = line.strip()
                    if line == '' or line.startswith('#'):
                        continue
                    try:
                        tag, value = [word.strip() for word in line.split('=')]
                    except ValueError:
                        print("")
                        print ('Invalid format:\n' + line)
                        raise UserWarning('Invalid model_file_index file {}.'.format(file_path))
                    if tag in self.parameters.keys():
                        raise UserWarning('Multiple entries for tag {} in {}.'.format(tag,file_path))

                    if (tag in self.obligatory_tags):
                        self.parameters[tag] = value
                    else:
                        print("")
                        print('WARNING: ignoring unknown parameter '
                              + tag + '.')
        except IOError:
            print("")
            print('Could not find model_file_index ' + file_path + '.')
            raise UserWarning('Invalid model_file_index file.')

        # check that all obligatory tags have been found
        missing_parameters = [tag for tag in self.obligatory_tags
                              if tag not in self.parameters.keys()]
        if len(missing_parameters) > 0:
            print("")
            print('Following tags are missing: '
                  + ', '.os.path.join(missing_parameters))
            raise UserWarning('Invalid model_file_index file {}.'.format(file_path))

    def write_to_file(self,file_path):
        with open(file_path, 'w') as f:
            for tag, value in self.parameters.items():
                f.write('{} = {}\n'.format(tag, value))
