"""Module importing user data needed to build the model."""

# python 2/3 compatibility
from __future__ import division, print_function, absolute_import

# global imports
import os.path

# local imports
from rba.prerba.pipeline_parameters import PipelineParameters
from rba.prerba.default_data import DefaultData
from rba.prerba import sbml_filter
from rba.prerba.uniprot_data import UniprotData
from rba.prerba.manual_annotation import (
    CuratedMetabolites, CuratedMacrocomponents
    )
from rba.prerba.protein_data import ProteinData
from rba.prerba.uniprot_importer import UniprotImporter
from rba.prerba.fasta_parser import FastaParser
from rba.prerba import protein_export


class UserData(object):

    def __init__(self, parameter_file):
        """Read data stored in filed described in parameters."""
        self._parameters = PipelineParameters(parameter_file).parameters
        self.default = DefaultData()
        print('Importing SBML data...')
        self.sbml_data = sbml_filter.SBMLFilter(
            self._sbml_file(), external_ids=self._external_ids()
            )
        print('Importing Uniprot data and manual annotation...')
        create_uniprot(self.input_path('uniprot.csv'), self._organism_id())
        self.protein_data = ProteinData(self._input_dir())
        self._initialize_gene_to_enzyme_mapping()
        self._initialize_average_protein()
        self.protein_data.update_helper_files()
        known_species = set([s.id for s in self.sbml_data.species])
        self.macrocomponents = CuratedMacrocomponents(
            self._input_dir(), known_species
            ).data
        self.metabolite_map = self.build_metabolite_map(
            self._input_dir(), known_species, self.cofactors()
            )
        self.rna_data = read_trnas(
            self.input_path('trnas.fasta'), self.metabolite_map
            )
        self.ribosome = FastaParser(self.input_path('ribosome.fasta')).entries
        self.chaperone = FastaParser(
            self.input_path('chaperones.fasta')
            ).entries

    def _sbml_file(self):
        return self.input_path(self._parameters['SBML_FILE'])

    def input_path(self, filename):
        return os.path.join(self._input_dir(), filename)

    def _input_dir(self):
        return self._parameters['INPUT_DIR']

    def _external_ids(self):
        line = self._parameters.get('EXTERNAL_COMPARTMENTS', None)
        if line is None:
            return []
        return [e.strip() for e in line.split(',')]

    def _organism_id(self):
        return self._parameters['ORGANISM_ID']

    def _initialize_gene_to_enzyme_mapping(self):
        self.enzymatic_proteins = {}
        self.protein_reference = {}
        for g in self._sbml_enzymatic_genes():
            protein, reference = self.protein_data.protein_and_reference(g)
            if protein:
                self.enzymatic_proteins[g] = protein
                self.protein_reference[g] = reference
            else:
                # spontaneous reaction or unknown protein
                if reference:
                    self.protein_reference[g] = reference

    def _sbml_enzymatic_genes(self):
        result = []
        for enzyme in self.sbml_data.enzymes:
            result += [g for g in enzyme if g != '']
        return list(set(result))

    def _initialize_average_protein(self):
        # we remove non-standard amino acids from the average composition
        average_protein = self.protein_data.average_composition()
        self.average_protein = {aa: sto for aa, sto in average_protein.items()
                                if aa in self.default.metabolites.aas}

    def output_dir(self):
        return self._parameters['OUTPUT_DIR']

    def export_proteins(self, filename):
        protein_export.export_proteins(
            self.input_path(filename), self.enzymatic_proteins
            )

    def sbml_species(self):
        return self.sbml_data.species

    def sbml_reactions(self):
        return self.sbml_data.reactions

    def enzyme_composition(self):
        return self.sbml_data.enzymes

    def has_membrane_enzyme(self, reaction):
        return self.sbml_data.has_membrane_enzyme[reaction]

    def imported_metabolites(self, reaction):
        return self.sbml_data.imported_metabolites.get(reaction, [])

    def external_metabolites(self):
        return self.sbml_data.external_metabolites

    def compartments(self):
        return self.protein_data.compartments()

    def compartment(self, id_):
        return self.protein_data.compartment(id_)

    def average_protein_id(self, compartment):
        """Return identifier of average protein in given compartment."""
        return self.protein_data.average_protein_id(compartment)

    def build_metabolite_map(self, input_dir, known_species, cofactors):
        """
        Map internal keys for metabolites with user-defined SBML ids.

        Parameters
        ----------
        input_dir : str
            Path to directory containing helper files.
        known_species : list of str
            List of valid metabolite identifiers.
        cofactors : list of rba.xml.uniprot_data.Cofactor
            List of cofactors.

        Returns
        -------
        dict
            Dictionary where keys are internal metabolite identifiers and
            values are corresponding user identifiers.

        """
        # metabolites to retrieve
        keys, names = self.default.metabolites.process_metabolites()
        cofactor_info = {}
        for cofactor in cofactors:
            cofactor_info.setdefault(cofactor.chebi, cofactor.name)
        keys += list(cofactor_info)
        names += list(cofactor_info.values())
        # look for metabolites
        curated_data = CuratedMetabolites(input_dir, known_species)
        sbml_lookup = {s.split('_', 1)[1].lower(): s for s in known_species}
        for key, name in zip(keys, names):
            # if curated data is available use it,
            # otherwise try to find sbml id using standard name
            if key not in curated_data.data:
                sbml_id = sbml_lookup.get((key + '_c').lower(), None)
                conc = self.default.metabolites.concentration.get(key, 0)
                curated_data.append(key, name, sbml_id, conc)
        curated_data.update_file()
        return curated_data.data

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
        for protein in self.enzymatic_proteins.values():
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


def ntp_composition(sequence):
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


def create_uniprot(input_file, organism_id):
    """
    Create uniprot file if necessary.

    Parameters
    ----------
    input_file : str
        File containing Uniprot data (or where it should be created).
    organism_id : int
        Uniprot organism id to fetch.

    """
    if not os.path.isfile(input_file):
        print('Could not find uniprot file. Downloading most recent'
              ' version...')
        raw_data = UniprotImporter(organism_id).data
        if len(raw_data) == 0:
            raise UserWarning(
                'Invalid organism, could not retrieve Uniprot data. '
                )
        with open(input_file, 'w') as f:
            f.write(raw_data)


def read_trnas(filename, metabolite_map):
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
        user_metabolite = metabolite_map.get(rna.id.upper(), None)
        if user_metabolite and user_metabolite.sbml_id:
            user_id = user_metabolite.sbml_id
        else:
            user_id = rna.id
        rna_list = sequence_list.setdefault(user_id, [])
        rna_list.append(rna.sequence)
    average_comp = {}
    for id_, seq in sequence_list.items():
        comp = ntp_composition(''.join(seq))
        average_comp[id_] = {k: v/len(seq) for k, v in comp.items()}
    return average_comp
