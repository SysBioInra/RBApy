"""Module importing user data needed to build the model."""

# python 2/3 compatibility
from __future__ import division, print_function, absolute_import

# global imports
import os.path
import sys

# local imports
from rba.prerba.pipeline_parameters import PipelineParameters
from rba.prerba.default_data import DefaultData
from rba.prerba import sbml_data
from rba.prerba.protein_data import ProteinData
from rba.prerba.uniprot_importer import create_uniprot_if_absent
from rba.prerba.manual_annotation import (
    CuratedMetabolites, CuratedMacrocomponents
)
from rba.prerba.user_machinery import UserMachinery
from rba.prerba.fasta_parser import RbaFastaParser
from rba.prerba import protein_export


class UserData(object):
    """Data contained in files provided by the user."""

    def __init__(self, parameter_file, verbose=False):
        """Read data stored in filed described in parameters."""
        self.verbose = verbose
        self._parameters = PipelineParameters(parameter_file).parameters
        self.default = DefaultData()
        self._import_sbml_data()
        self._import_uniprot_data()
        self._import_manual_annotation()

    def _import_sbml_data(self):
        if self.verbose:
            print('  Importing SBML data ...', end='')
            sys.stdout.flush()
        self.sbml_data = sbml_data.SbmlData(
            self.input_path(self._parameters['SBML_FILE']),
            external_ids=self._external_ids(),
            interface_id=self._interface_ids()
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
        create_uniprot_if_absent(self.input_path('uniprot.csv'),
                                 self._organism_id())
        self.protein_data = ProteinData(self._input_dir())
        self._retrieve_enzymatic_proteins()
        self.protein_data.update_helper_files()
        if self.verbose:
            print(' done')

    def _organism_id(self):
        return self._parameters['ORGANISM_ID']

    def _retrieve_enzymatic_proteins(self):
        self.enzymatic_proteins = []
        for g in self._sbml_enzymatic_genes():
            protein = self.protein_data.create_protein_from_gene_id(g)
            if protein:
                self.enzymatic_proteins.append(protein)

    def _sbml_enzymatic_genes(self):
        result = []
        for enzyme in self.sbml_data.enzymes:
            result += [g for g in enzyme.gene_association if g != '']
        return list(set(result))

    def _import_manual_annotation(self):
        if self.verbose:
            print('  Importing manual annotation ...', end='')
            sys.stdout.flush()
        known_species = self._sbml_species_ids()
        self.macrocomponents = CuratedMacrocomponents(
            self._input_dir(), known_species
        ).data
        self.metabolite_map = self._build_metabolite_map()
        self.trnas = self._read_trnas(self.input_path('trnas.fasta'))
        self.ribosome = UserMachinery(self.input_path('ribosome.fasta'),
                                      self.protein_data)
        self.chaperone = UserMachinery(self.input_path('chaperones.fasta'),
                                       self.protein_data)
        if self.verbose:
            print(' done')

    def _sbml_species_ids(self):
        return set([s.id for s in self.sbml_data.species])

    def _build_metabolite_map(self):
        """Map internal keys for metabolites with user-defined SBML ids."""
        known_species = self._sbml_species_ids()
        curated_data = CuratedMetabolites(self._input_dir(), known_species)
        sbml_lookup = {s.split('_', 1)[1].lower(): s for s in known_species}
        for id_, name in zip(*self._internal_species_ids_and_names()):
            if id_ not in curated_data.data:
                # id_ not mapped in curation file: add new entry
                sbml_id = sbml_lookup.get((id_ + '_c').lower(), None)
                conc = self.default.metabolites.concentration.get(id_, 0)
                curated_data.append(id_, name, sbml_id, conc)
        curated_data.update_file()
        return curated_data.data

    def _internal_species_ids_and_names(self):
        keys, names = self.default.metabolites.process_metabolites()
        cofactor_info = {}
        for cofactor in self.cofactors():
            cofactor_info.setdefault(cofactor.chebi, cofactor.name)
        keys += list(cofactor_info)
        names += list(cofactor_info.values())
        return keys, names

    def _read_trnas(self, filename):
        """Read trnas in fasta file."""
        trna_data = RbaFastaParser(filename).rnas
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

    def average_protein_id(self, compartment):
        """Return identifier of average protein in given compartment."""
        return self.protein_data.average_protein_id(compartment)

    def transport_enzymes(self):
        return (e for e in self.sbml_data.enzymes if e.is_transporter)

    def metabolite_targets(self):
        result = [(id_, m.sbml_id, m.concentration)
                  for id_, m in self.metabolite_map.items()
                  if m.sbml_id and m.concentration]
        result += [(id_, id_, conc)
                   for id_, conc in self.macrocomponents.items()]
        return result

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

    def sbml_enzymes(self):
        result = self.sbml_data.enzymes
        for enzyme in result:
            enzyme.composition = self._build_enzyme_composition(
                enzyme.gene_association
            )
        return result

    def _build_enzyme_composition(self, composition):
        result = []
        for gene in composition:
            reference = self.protein_data.reference(gene)
            if reference:
                result.append(reference)
        return result

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
        for protein in self.enzymatic_proteins:
            for c in protein.cofactors:
                if c.chebi not in known_ids:
                    cofactors.append(c)
                    known_ids.add(c.chebi)
        return cofactors
