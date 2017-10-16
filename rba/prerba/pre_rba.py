"""Module defining PreRba class."""

# python 2/3 compatibility
from __future__ import division, print_function, absolute_import

# global imports
from collections import namedtuple
import itertools
import os.path
import pandas

# data processing imports
from rba.prerba.pipeline_parameters import PipelineParameters
from rba.prerba import sbml_filter
from rba.prerba.uniprot_data import UniprotData
from rba.prerba.manual_annotation import ManualAnnotation
from rba.prerba.protein_data import ProteinData
from rba.prerba.uniprot_importer import UniprotImporter
from rba.prerba.fasta_parser import FastaParser

# default data imports
from rba.prerba.default_data import DefaultData
from rba.prerba import default_processes

# model imports
import rba.xml
from rba.rba_model import RbaModel

Metabolite = namedtuple('Metabolite', 'name sbml_id concentration')


class PreRba(object):
    """
    Class generating RBA model from SBML file and helper files.

    Attributes
    ----------
    model: rba.rba_model.RbaModel
        RBA model generated.
    metabolite_map: dict
        Map fromt internal metabolite names to sbml identifiers.

    """

    def __init__(self, parameter_file):
        """
        Build from parameter_file.

        Parameters
        ----------
        parameter_file: str
            standard pipeline input file (containing e.g.
            sbml location, organism id).

        """
        parameters = PipelineParameters(parameter_file).parameters
        print('Importing default data...')
        self.default = DefaultData()
        print('Importing user data...')
        self._read_data(parameters)
        print('Building model...')
        self.model = self.build_model()
        self.model.output_dir = parameters['OUTPUT_DIR']
        print('Done.')

    def _read_data(self, parameters):
        """Read data stored in filed described in parameters."""
        organism_id = parameters['ORGANISM_ID']
        sbml_file = parameters['SBML_FILE']
        input_dir = parameters['INPUT_DIR']
        external_line = parameters.get('EXTERNAL_COMPARTMENTS', None)
        if external_line is not None:
            external_ids = external_line.split(',')
            external_ids = [e.strip() for e in external_compartments]
        else:
            external_ids = []
        print('Importing SBML data...')
        self.sbml_data = sbml_filter.SBMLFilter(
            os.path.join(input_dir, sbml_file), external_ids=external_ids
            )
        genes_to_retrieve = []
        for enzyme in self.sbml_data.enzymes:
            genes_to_retrieve += [g for g in enzyme]
        genes_to_retrieve = list(set(genes_to_retrieve))
        print('Importing Uniprot data and manual annotation...')
        uniprot = open_uniprot(os.path.join(input_dir, 'uniprot.csv'),
                               organism_id)
        manual = ManualAnnotation(input_dir)
        self.protein_data, self.protein_stoichiometry, self.location_map \
            = protein_data(uniprot, manual, genes_to_retrieve)
        average_protein = uniprot.average_protein_composition()
        # remove non-standard amino acids
        self.average_protein = {aa: sto for aa, sto in average_protein.items()
                                if aa in self.default.metabolites.aas}
        self.rna_data = read_trnas(
            os.path.join(input_dir, 'trnas.fasta'),
            )
        self.ribosome = FastaParser(
            os.path.join(input_dir, 'ribosome.fasta')
            ).entries
        self.chaperone = FastaParser(
            os.path.join(input_dir, 'chaperones.fasta')
            ).entries
        known_species = set([s.id for s in self.sbml_data.species])
        self.metabolite_map = self.build_metabolite_map(
            known_species, manual.metabolites, self.cofactors()
            )
        self.macrocomponents = {}
        for met, flux in manual.macrocomponents.data.data.values:
            if met in known_species:
                self.macrocomponents[met] = float(flux)
            else:
                print('WARNING: in file {}: {} is not a valid SBML species. '
                      'Line will be ignored.'
                      .format(manual.macrocomponents.filename, met))

    def build_model(self):
        """
        Build RBA model according to data read.

        Returns
        -------
        rba.rba_model.RbaModel
            Full RBA model.

        """
        model = RbaModel()
        model.metabolism = self.build_metabolism()
        model.parameters = self.build_parameters()
        model.proteins = self.build_proteins()
        model.rnas = self.build_rnas()
        model.dna = self.build_dna()
        model.processes = self.build_processes()
        model.enzymes = self.build_enzymes()
        concentration = self.default.activity.medium_concentration
        for metab in self.sbml_data.external_metabolites:
            # !!! we identify metabolites by their prefix !!!
            # the idea is that M_glc_p and M_glc_e will be seen
            # as the same metabolite M_glc
            model.medium[metab.rsplit('_', 1)[0]] = concentration
        return model

    def build_metabolism(self):
        """
        Build metabolism part of RBA model.

        Returns
        -------
        rba.xml.RbaMetabolism
            RBA metabolism model in XML format.

        """
        metabolism = rba.xml.RbaMetabolism()
        metabolism.species = self.sbml_data.species
        metabolism.reactions = self.sbml_data.reactions
        # add atpm reaction
        reaction = rba.xml.Reaction(self.default.atpm_reaction, False)
        for m in ['ATP', 'H2O']:
            id_ = self.metabolite_map[m].sbml_id
            if id_:
                reaction.reactants.append(rba.xml.SpeciesReference(id_, 1))
        for m in ['ADP', 'H', 'Pi']:
            id_ = self.metabolite_map[m].sbml_id
            if id_:
                reaction.products.append(rba.xml.SpeciesReference(id_, 1))
        metabolism.reactions.append(reaction)
        for c in self.compartments():
            metabolism.compartments.append(rba.xml.Compartment(c))
        return metabolism

    def compartments(self):
        """
        Return list of compartment identifiers.

        Returns
        -------
        list
            List of compartment identifiers.

        """
        compartment_ids = []
        for uniprot_id, user_id in self.location_map.items():
            if not pandas.isnull(user_id):
                compartment_ids.append(user_id)
            else:
                compartment_ids.append(uniprot_id.replace(' ', '_'))
        return list(set(compartment_ids))

    def compartment(self, uniprot_id):
        """
        Return compartment corresponding to a given uniprot identifier.

        Parameters
        ----------
        uniprot_id: str
            Valid uniprot compartment identifier (e.g. Cytoplasm, Secreted).

        Returns
        -------
        str
            Model identifier associated with compartment.

        """
        user_id = self.location_map[uniprot_id]
        if pandas.isnull(user_id):
            return uniprot_id.replace(' ', '_')
        else:
            return user_id

    def build_parameters(self):
        """
        Build parameter part of RBA model.

        Returns
        -------
        rba.xml.RbaParameters
            RBA parameter model in XML format.

        """
        parameters = rba.xml.RbaParameters()
        def_params = self.default.parameters
        # target densities
        cytoplasm = self.compartment('Cytoplasm')
        external = self.compartment('Secreted')
        other_compartments = self.compartments()
        other_compartments.remove(cytoplasm)
        other_compartments.remove(external)
        constraints = parameters.target_densities
        aggregates = parameters.aggregates
        # always start with cytoplasm
        c_density = rba.xml.TargetDensity(cytoplasm)
        c_density.upper_bound = def_params.cytoplasm_density
        constraints.append(c_density)
        # now treat the rest
        for user_id in other_compartments:
            # create aggregated function for density
            id_ = user_id + '_density'
            new_aggregate = rba.xml.Aggregate(id_, 'multiplication')
            for ref in ['amino_acid_concentration',
                        def_params.protein_fraction_id(user_id)]:
                new_aggregate.function_references.append(
                    rba.xml.FunctionReference(ref)
                    )
            aggregates.append(new_aggregate)
            # create constraint
            new_density = rba.xml.TargetDensity(user_id)
            new_density.upper_bound = id_
            constraints.append(new_density)

        # functions
        fns = parameters.functions
        # protein fractions
        fns.append(def_params.protein_fraction(
            cytoplasm, def_params.cytoplasm_fraction
            ))
        fns.append(def_params.protein_fraction(
            external, def_params.secreted_fraction
            ))
        other_fraction = def_params.other_fraction / len(other_compartments)
        for cpt in other_compartments:
            fns.append(def_params.protein_fraction(cpt, other_fraction))
        # non enzymatic fractions
        fns.append(def_params.non_enzymatic_fraction_cytoplasm(cytoplasm))
        fns.append(def_params.non_enzymatic_fraction_secreted(external))
        for comp in other_compartments:
            fns.append(def_params.non_enzymatic_fraction_other(comp))
        # protein length
        fns.append(def_params.inverse_average_protein_length(
            sum(self.average_protein.values())
            ))
        # add all other functions
        for function in def_params.functions:
            fns.append(function)
        # add aggregates
        for agg in default_processes.aggregates(def_params,
                                                self.compartments()):
            aggregates.append(agg)
        return parameters

    def cofactors(self):
        """
        Extract protein cofactors.

        Returns
        -------
        list
            List of cofactors associated with protein in the model. Each
            cofactor is represented only once.

        """
        cofactors = []
        known_ids = set()
        for protein in self.protein_data.values():
            for c in protein.cofactors:
                if c.chebi not in known_ids:
                    cofactors.append(c)
                    known_ids.add(c.chebi)
        return cofactors

    def aa_composition(self, sequence):
        """
        Translate sequence into amino acid composition.

        Parameters
        ----------
        sequence : str
            Protein sequence in one-letter format.

        Returns
        -------
        dict
            Dictionary where keys are amino acid identifiers and values the
            number of times the amino acid appeared in the sequence.

        """
        return composition(sequence, self.default.metabolites.aas)

    def build_proteins(self):
        """
        Build protein part of RBA model.

        Returns
        -------
        rba.xml.RbaMacromolecules
            RBA protein model in XML format.

        """
        proteins = rba.xml.RbaMacromolecules()
        # components
        for aa in self.default.metabolites.aas:
            proteins.components.append(
                rba.xml.Component(aa, '', 'amino_acid', 1)
                )
        for c in self.cofactors():
            proteins.components.append(
                rba.xml.Component(c.chebi, c.name, 'cofactor', 0)
                )
        # user proteins
        for gene_name, protein in self.protein_data.items():
            comp = self.aa_composition(protein.sequence)
            for cofactor in protein.cofactors:
                comp[cofactor.chebi] = cofactor.stoichiometry
            proteins.macromolecules.append(
                rba.xml.Macromolecule(gene_name, protein.location, comp)
                )
        # average proteins
        for compartment in self.compartments():
            proteins.macromolecules.append(
                rba.xml.Macromolecule('average_protein_' + compartment,
                                      compartment, self.average_protein)
                )
        # machinery proteins
        cytoplasm = self.compartment('Cytoplasm')
        for molecule in itertools.chain(self.ribosome, self.chaperone):
            if molecule.set_name == 'protein':
                proteins.macromolecules.append(
                    rba.xml.Macromolecule(
                        molecule.id, cytoplasm,
                        self.aa_composition(molecule.sequence)
                    )
                )
        return proteins

    def ntp_composition(self, sequence):
        """
        Translate sequence into ntp composition.

        Parameters
        ----------
        sequence : str
            DNA or RNA sequence. If a DNA sequence is provided, all T will
            be converted to U.

        Returns
        -------
        dict
            Dictionary where keys are ntp identifiers and values the
            number of times the ntps appeared in the sequence.

        """
        comp = composition(sequence, 'ACGTU')
        comp['U'] += comp.pop('T')
        return comp

    def build_rnas(self):
        """
        Build RNA part of RBA model.

        Returns
        -------
        rba.xml.RbaMacromolecules
            RBA RNA model in XML format.

        """
        rnas = rba.xml.RbaMacromolecules()
        # components
        rnas.components.append(
            rba.xml.Component('A', 'Adenosine residue', 'Nucleotide', 2.9036)
            )
        rnas.components.append(
            rba.xml.Component('C', 'Cytosine residue', 'Nucleotide', 2.7017)
            )
        rnas.components.append(
            rba.xml.Component('G', 'Guanine residue', 'Nucleotide', 3.0382)
            )
        rnas.components.append(
            rba.xml.Component('U', 'Uramine residue', 'Nucleotide', 2.7102)
            )
        # user rnas
        cytoplasm = self.compartment('Cytoplasm')
        for rna_id, composition in self.rna_data.items():
            user_rna = self.metabolite_map.get(rna_id.upper(), None)
            if user_rna and user_rna.sbml_id:
                rnas.macromolecules.append(rba.xml.Macromolecule(
                    user_rna.sbml_id, cytoplasm, composition
                    ))
        # average RNA
        rnas.macromolecules.append(rba.xml.Macromolecule(
            self.default.metabolites.mrna, cytoplasm,
            {'A': 0.2818, 'C': 0.2181, 'G': 0.2171, 'U': 0.283}
            ))
        # machinery rnas
        for molecule in itertools.chain(self.ribosome, self.chaperone):
            if molecule.set_name == 'rna':
                rnas.macromolecules.append(rba.xml.Macromolecule(
                    molecule.id, cytoplasm,
                    self.ntp_composition(molecule.sequence)
                    ))
        return rnas

    def build_dna(self):
        """
        Build DNA part of RBA model.

        Returns
        -------
        rba.xml.RbaMacromolecules
            RBA DNA model in XML format.

        """
        dna = rba.xml.RbaMacromolecules()
        # components
        comp_a = rba.xml.Component('A', 'Adenosine residue', 'Nucleotide', 0)
        comp_c = rba.xml.Component('C', 'Cytosine residue', 'Nucleotide', 0)
        comp_g = rba.xml.Component('G', 'Guanine residue', 'Nucleotide', 0)
        comp_t = rba.xml.Component('T', 'Thymine residue', 'Nucleotide', 0)
        for comp in [comp_a, comp_c, comp_g, comp_t]:
            dna.components.append(comp)
        # average DNA
        dna.macromolecules.append(
            rba.xml.Macromolecule(
                self.default.metabolites.dna, self.compartment('Cytoplasm'),
                {'A': 0.2818, 'C': 0.2181, 'G': 0.2171, 'T': 0.283}
                )
            )
        return dna

    def build_enzymes(self):
        """
        Build enzyme part of RBA model.

        Returns
        -------
        rba.xml.RbaEnzymes
            RBA enzyme model in XML format.

        """
        enzymes = rba.xml.RbaEnzymes()
        cytoplasm = self.compartment('Cytoplasm')
        # efficiency functions
        # for the moment simply put default value
        default_function = rba.xml.Function('default', 'constant', {})
        enzymes.efficiency_functions.append(default_function)

        # user enzymes
        default_activity = self.default.activity
        def_param = {'CONSTANT': default_activity.catalytic_activity}
        def_efficiency = rba.xml.EnzymeEfficiency('default', def_param)
        tra_param = {'CONSTANT': default_activity.transporter_activity}
        tra_efficiency = rba.xml.EnzymeEfficiency('default', tra_param)
        mm_parameters = {'Km': default_activity.import_Km,
                         'kmax': default_activity.import_kmax}
        reaction_ids = [r.id for r in self.sbml_data.reactions]
        for reaction, comp in zip(reaction_ids, self.sbml_data.enzymes):
            # base information
            enzyme = rba.xml.Enzyme(reaction)
            enzymes.enzymes.append(enzyme)
            # machinery composition
            reactants = enzyme.machinery_composition.reactants
            for gene in comp:
                name_sto = self.protein_stoichiometry.get(gene, None)
                if name_sto:
                    reactants.append(rba.xml.SpeciesReference(*name_sto))
            # base enzymatic activity
            efficiencies = enzyme.enzymatic_activity.enzyme_efficiencies
            if self.sbml_data.has_membrane_enzyme[reaction]:
                efficiencies.append(tra_efficiency)
            else:
                efficiencies.append(def_efficiency)
            # transport activity
            transport = enzyme.enzymatic_activity.transporter_efficiency
            imported = self.sbml_data.imported_metabolites.get(reaction, [])
            for m in imported:
                transport.append(
                    rba.xml.Function('', 'michaelisMenten', mm_parameters, m)
                    )

        # atpm enzyme
        enzyme = rba.xml.Enzyme(self.default.atpm_reaction)
        efficiency = rba.xml.EnzymeEfficiency(
            'default', {'CONSTANT': self.default.activity.catalytic_activity}
            )
        enzyme.enzymatic_activity.enzyme_efficiencies.append(efficiency)
        enzymes.enzymes.append(enzyme)
        return enzymes

    def build_processes(self):
        """
        Build process part of RBA model.

        Returns
        -------
        rba.xml.RbaProcesses
            RBA process model in XML format.

        """
        processes = rba.xml.RbaProcesses()
        def_proc = default_processes.DefaultProcesses(self.default,
                                                      self.metabolite_map)
        compartments = self.compartments()
        # processes
        proc_list = processes.processes
        proc_list.append(def_proc.translation(
            {m.id: m.stoichiometry for m in self.ribosome}, compartments
            ))
        proc_list.append(def_proc.folding(
            {m.id: m.stoichiometry for m in self.chaperone}
            ))
        proc_list.append(def_proc.transcription())
        proc_list.append(def_proc.replication())
        proc_list.append(def_proc.rna_degradation())
        proc_list.append(def_proc.metabolite_production())
        proc_list.append(def_proc.macrocomponents(self.macrocomponents))
        proc_list.append(def_proc.maintenance_atp(self.default.atpm_reaction))
        # component maps
        map_list = processes.component_maps
        map_list.append(def_proc.translation_map(self.cofactors()))
        map_list.append(def_proc.folding_map())
        map_list.append(def_proc.transcription_map())
        map_list.append(def_proc.rna_degradation_map())
        map_list.append(def_proc.replication_map())
        return processes

    def build_metabolite_map(self, known_species, curated_data, cofactors):
        """
        Map internal keys for metabolites with user-defined SBML ids.

        Parameters
        ----------
        known_species : list of str
            List of valid metabolite identifiers.
        curated_data : rba.xml.manual_annotation.CurationFile
            Manual annotation file containing metabolite information.
        cofactors : list of rba.xml.uniprot_data.Cofactor
            List of cofactors.

        Returns
        -------
        dict
            Dictionary where keys are internal metabolite identifiers and
            values are corresponding user identifiers.

        """
        # build dictionary with curated metabolites
        curated_metabolites = {}
        for id_, name, sbml_id, conc in curated_data.data.data.values:
            if pandas.isnull(sbml_id) or sbml_id in known_species:
                if pandas.isnull(sbml_id):
                    sbml_id = None
                if pandas.isnull(conc) or conc == '':
                    conc = 0
                curated_metabolites[id_] = Metabolite(
                    name, sbml_id, float(conc)
                    )
            else:
                print('ERROR: {} is not a valid metabolite id.'
                      .format(sbml_id))

        # retrieve items
        metabolites = {}
        data_to_cure = []
        keys, names = self.default.metabolites.process_metabolites()
        cofactor_info = {}
        for cofactor in cofactors:
            cofactor_info.setdefault(cofactor.chebi, cofactor.name)
        keys += list(cofactor_info)
        names += list(cofactor_info.values())

        sbml_lookup = {s.split('_', 1)[1].lower(): s for s in known_species}
        for key, name in zip(keys, names):
            # if curated data is available use it,
            # otherwise try to find sbml id using standard name
            met = curated_metabolites.get(key, None)
            if met:
                metabolites[key] = met
            else:
                sbml_id = sbml_lookup.get((key + '_c').lower(), None)
                conc = self.default.metabolites.concentration.get(key, 0)
                data_to_cure.append((key, name, sbml_id, conc))
                metabolites[key] = Metabolite(name, sbml_id, conc)

        # store data to cure
        if data_to_cure:
            curated_data.data.add(data_to_cure)
            curated_data.data.write(curated_data.filename)
        return metabolites


def open_uniprot(input_file, organism_id):
    """
    Open uniprot file (create it if necessary).

    Parameters
    ----------
    input_file : str
        File containing Uniprot data (or where it should be created).
    organism_id : int
        Uniprot organism id to fetch.

    Returns
    -------
    rba.xml.uniprot_data.UniprotData
        Uniprot data.

    """
    if not os.path.isfile(input_file):
        print('Could not find uniprot file. Downloading most recent'
              ' version...')
        raw_data = UniprotImporter(organism_id).data
        if len(raw_data) == 0:
            raise UserWarning(
                'Invalid organism, could not retrieve Uniprot data. '
                'Interrupting execution.'
                )
        with open(input_file, 'w') as f:
            f.write(raw_data)
    return UniprotData(input_file)


def protein_data(uniprot, manual, sbml_ids):
    """
    Retrieve protein information.

    Parameters
    ----------
    uniprot : rba.xml.uniprot_data.UniprotData
        Uniprot data.
    manual : rba.xml.manual_annotation.ManualAnnotation
        Manual annotation.
    sbml_ids : list of str
        List of genes to retrieve

    Returns
    -------
    dict
        Base protein info. Keys are genes id provided as input, values are
        rba.xml.protein_data.Protein structures.
    dict
        Stoichiometry info. Keys are genes id provided as input, values are
        tuples (protein_id, stoichiometry).
    dict
        Location map. Keys are uniprot ids for locations, values are
        corresponding user locations.

    """
    sbml_ids = [s for s in sbml_ids if s != '']
    sbml_to_user = {}
    for sbml_id, user_id in manual.unknown_proteins.data.data.values:
        sbml_to_user[sbml_id] = sbml_id if pandas.isnull(user_id) else user_id
    protein_data = ProteinData(uniprot, manual)
    prot_info = {}
    stoichiometry = {}
    not_found = []
    invalid_user_ids = []
    average_protein = 'average_protein_' + protein_data.default_location
    for sbml_id in sbml_ids:
        user_id = sbml_to_user.get(sbml_id, 0)
        uniprot_id = None
        if user_id == 0:
            uniprot_id = uniprot.entry(sbml_id)
            if not uniprot_id:
                stoichiometry[sbml_id] = (average_protein, 1)
                not_found.append((sbml_id, average_protein))
        else:
            # check whether id was mapped to a spontaneous reaction,
            # an average protein or a new gene id
            if user_id == '' or pandas.isnull(user_id):
                pass
            elif user_id.startswith('average_protein_'):
                stoichiometry[sbml_id] = (user_id, 1)
            else:
                uniprot_id = uniprot.entry(user_id)
                if not uniprot_id:
                    invalid_user_ids.append(user_id)
        if uniprot_id:
            prot = protein_data.find_protein(uniprot_id)
            prot_info[sbml_id] = prot
            stoichiometry[sbml_id] = (sbml_id, prot.stoichiometry)
    protein_data.update_manual_annotation(manual)
    if not_found:
        manual.unknown_proteins.data.add(not_found)
        manual.unknown_proteins.data.write(manual.unknown_proteins.filename)
    if invalid_user_ids:
        print('WARNING: {} are invalid protein identifiers. '
              'Check data provided provided in {}.'
              .format(', '.join(invalid_user_ids),
                      manual.unknown_proteins.filename))
    return prot_info, stoichiometry, protein_data.location_map


def composition(sequence, alphabet):
    """
    Compute composition of sequence with given alphabet.

    Parameters
    ----------
    sequence : str
        Sequence to decompose.
    alphabet : str
        Letters to count in the original sequence. A letter of the
        sequence that is not in alphabet will be ignored.

    Returns
    -------
    dict
        Keys are letters of the alphabet and values the corresponding counts
        in the sequence.

    """
    comp = dict.fromkeys(alphabet, 0)
    for n in sequence:
        try:
            comp[n] += 1
        except KeyError:
            pass
    return comp


def read_trnas(filename):
    """
    Read trnas in fasta file.

    Parameters
    ----------
    filename : str
        File containing RNA information.

    Returns
    -------
    dict
        Keys are identifiers and values are the composition of RNAs.
        A composition is a dictionary where keys are NTPs and values are the
        average counts for the NTP over RNAs sharing the same identifier.

    """
    # read all real trnas (as described in fasta files)
    trna_data = FastaParser(filename).entries
    # map real trnas to user trnas
    # for example, user may agregate all trnas into a single metabolite
    # in this case, we take an average composition for a trna
    sequence_list = {}
    for rna in trna_data:
        rna_list = sequence_list.setdefault(rna.id.upper(), [])
        rna_list.append(rna.sequence)
    average_comp = {}
    for id_, seq in sequence_list.items():
        comp = composition(''.join(seq), 'ACGTU')
        comp['U'] += comp.pop('T')
        average_comp[id_] = {k: v/len(seq) for k, v in comp.items()}
    return average_comp