"""
Module defining Metabolite and MetaboliteData classes.
"""

# python 2/3 compatibility
from __future__ import division, print_function

# global imports
import os.path
from collections import namedtuple
import pandas

# local imports
from rba.prerba.curation_data import CurationData

# Class used to store metabolite information
Metabolite = namedtuple('Metabolite', 'name sbml_id concentration')

class MetaboliteData(object):
    """
    Class used to parse/store metabolite information.

    Attributes:
        metabolites: dict matching internal metabolite keys, 
            with existing SBML identifiers.
    """

    _helper_file = 'metabolites.tsv'

    def __init__(self, known_species, default_metabolites,
                 cofactors, input_dir='.'):
        """
        Constructor.

        Args:
            known_species: identifiers of species defined in SBML file.
            default_metabolites: default metabolites as defined in 
                default_data.DefaultData.
            cofactors: list of known Cofactors.
            input_dir: path to folder containing metabolite files.
        """
        print('\nParsing metabolite information'
              '\n------------------------------')

        # read curation data
        curation_file = os.path.join(input_dir, self._helper_file)
        curation_data = CurationData(['ID', 'NAME', 'SBML ID', 'CONCENTRATION'])
        try:
            curation_data.read(curation_file)
            print('Helper file found.')
        except IOError:
            print('Helper file not found.')
        curated_metabolites = {}
        for id_, name, sbml_id, conc in curation_data.data.values:
            if pandas.isnull(sbml_id) or sbml_id in known_species:
                if pandas.isnull(sbml_id):
                    sbml_id = None
                if pandas.isnull(conc) or conc=='':
                    conc = 0
                curated_metabolites[id_] = Metabolite(name, sbml_id,
                                                      float(conc))
            else:
                print('ERROR: {} is not a valid metabolite id.'.format(sbml_id))

        # extract metabolite information
        self.metabolites, data_to_cure \
            = self._find_mapping(known_species, default_metabolites,
                                 cofactors, curated_metabolites)
        if data_to_cure:
            curation_data.add(data_to_cure)
            curation_data.write(curation_file)

        # raise warning if data was missing
        if curation_data.has_missing_information():
            print('WARNING: Several key metabolites could not be identified. '
                  'Please read file {}, check data and specify all missing '
                  'information. Unknown metabolites will be removed '
                  'from production targets for this run.'
                  .format(curation_file))

    @staticmethod
    def _ids_and_names(default_metabolites, cofactors):
        """
        Return internal ids and names of metabolites to retrieve.
        """
        # metabolites to retrieve
        keys = [] # internal id of metabolite
        names = [] # name of metabolite
        # methionine
        keys.append('MET')
        names.append('Methionine')
        # charged + uncharged trnas
        for aa in default_metabolites.aas:
            keys.append(default_metabolites.charged_trna_key(aa))
            names.append(default_metabolites.charged_trna_name(aa))
            keys.append(default_metabolites.uncharged_trna_key(aa))
            names.append(default_metabolites.uncharged_trna_name(aa))
        keys.append(default_metabolites.charged_trna_key
                    (default_metabolites.aa_fM))
        names.append(default_metabolites.charged_trna_name
                     (default_metabolites.aa_fM))
        # nucleotides
        for nucleotide in default_metabolites.nucleotides:
            keys.append(default_metabolites.ntp_key(nucleotide))
            names.append(default_metabolites.ntp_key(nucleotide))
            keys.append(default_metabolites.ndp_key(nucleotide))
            names.append(default_metabolites.ndp_key(nucleotide))
            keys.append(default_metabolites.nmp_key(nucleotide))
            names.append(default_metabolites.nmp_key(nucleotide))
        for nucleotide in default_metabolites.d_nucleotides:
            keys.append(default_metabolites.dntp_key(nucleotide))
            names.append(default_metabolites.dntp_key(nucleotide))
        # key metabolites
        for met_id, name in default_metabolites.key_metabolites.items():
            keys.append(met_id)
            names.append(name)
        # cofactors
        cofactor_info = {}
        for cofactor in cofactors:
            cofactor_info.setdefault(cofactor.id, cofactor.name)
        return keys + list(cofactor_info), names + list(cofactor_info.values())

    def _find_mapping(self, known_species, default_metabolites,
                      cofactors, curated_metabolites):
        """
        Map internal keys for metabolites with user-defined SBML ids.
        """
        # extract metabolite prefix
        met_prefix = known_species[0].split('_', 1)[0] + '_'
        assert all((m.startswith(met_prefix) for m in known_species))

        # retrieve items
        metabolites = {}
        data_to_cure = []
        keys, names = self._ids_and_names(default_metabolites, cofactors)
        sbml_lookup = {s.split('_', 1)[1].lower(): s for s in known_species}
        for key, name in zip(keys, names):
            # if curated data is available use it,
            # otherwise try to find sbml id using standard name
            met = curated_metabolites.get(key, None)
            if met:
                metabolites[key] = met
            else:
                # try simple standard name as sbml id
                sbml_id = sbml_lookup.get((key + '_c').lower(), None)
                conc = default_metabolites.concentration.get(key, 0)
                data_to_cure.append((key, name, sbml_id, conc))
                metabolites[key] = Metabolite(name, sbml_id, conc)
        return metabolites, data_to_cure
