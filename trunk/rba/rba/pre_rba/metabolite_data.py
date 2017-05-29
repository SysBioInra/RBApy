
import os.path
import re
import pandas

from . import rba_data

class Metabolite(object):
    def __init__(self, name, sbml_id, concentration):
        self.name = name
        self.sbml_id = sbml_id
        self.concentration = concentration

class MetaboliteData(object):
    """
    Class used to parse/store metabolite information.
    """
    def __init__(self, known_species, cofactors, input_dir = '.'):
        """
        Constructor.
        :param input_dir: path to folder containing metabolite files.
        :type input_dir: string
        """
        # constant values
        self._missing_string = ''
        self._missing_tag = '[MISSING]'
        self._data_file = os.path.join(input_dir, 'metabolites.tsv')
        
        # read curated data
        self._missing_information = False
        self.metabolites = {}
        self._read_curated_data(known_species)

        # extract metabolite information
        self._extract_metabolite_information(known_species, cofactors)
        
        # raise warning if data was missing
        if self._missing_information:
            print('\nWARNING: Several key metabolites could not be identified. '
                  'Please read file ' + self._data_file
                  + ', check data and specify all information tagged as '
                  + self._missing_tag + '. Unknown metabolites will be removed '
                  'from production targets for this run. Make sure to update '
                  'helper file for next execution.\n')

    def _read_curated_data(self, known_species):
        """
        Read file containing hand-curated information.
        """
        try:
            with open(self._data_file, 'rU') as input_stream:
                print('Found file with metabolite data. This file will be used '
                      'to identify key metabolites...')
                # skip header
                next(input_stream)
                # read lines
                for line in input_stream:
                    line = line.rstrip('\n')
                    [id_, name, sbml_id, conc] = line.split('\t')
                    if sbml_id == self._missing_tag: 
                        self._missing_information = True
                        sbml_id = self._missing_string
                    elif sbml_id not in known_species:
                        print('ERROR: ' + sbml_id +
                              ' is not a valid metabolite id.')
                        self._missing_information = True
                        sbml_id = self._missing_string
                    try:
                        conc = float(conc)
                    except ValueError:
                        conc = 0
                    self.metabolites[id_] = Metabolite(name, sbml_id, conc)
        except IOError:
            print('Could not find file with metabolite data...')

    def _write_curated_data(self):
        """
        Write file containing hand-curated metabolite information.
        """
        with open(self._data_file, 'w') as output_stream:
            output_stream.write('\t'.join(['ID', 'NAME', 'SBML ID',
                                           'CONCENTRATION']) + '\n')
            for id_, metab in self.metabolites.items():
                if metab.sbml_id == self._missing_string:
                    sbml_id = self._missing_tag
                else:
                    sbml_id = metab.sbml_id
                if metab.concentration != 0:
                    conc = str(metab.concentration)
                else:
                    conc = ''
                output_stream.write('\t'.join([id_, metab.name,
                                               sbml_id, conc]) + '\n')
        
    def _extract_metabolite_information(self, known_species, cofactors):
        curated_data_added = False

        # extract metabolite prefix
        met_prefix = known_species[0].split('_',1)[0] + '_'
        assert(all((m.startswith(met_prefix) for m in known_species)))

        # metabolites to retrieve
        keys = [] # internal id of metabolite
        names = [] # name of metabolite
        # methionine
        keys.append('MET')
        names.append('Methionine')
        # charged + uncharged trnas
        for aa in rba_data.aas:
            keys.append(rba_data.charged_trna_key(aa))
            names.append(rba_data.charged_trna_name(aa))
            keys.append(rba_data.uncharged_trna_key(aa))
            names.append(rba_data.uncharged_trna_name(aa))
        keys.append(rba_data.charged_trna_key(rba_data.aa_fM))
        names.append(rba_data.charged_trna_name(rba_data.aa_fM))
        # nucleotides
        for n in rba_data.nucleotides:
            keys.append(rba_data.ntp_key(n))
            names.append(rba_data.ntp_key(n))
            keys.append(rba_data.ndp_key(n))
            names.append(rba_data.ndp_key(n))
            keys.append(rba_data.nmp_key(n))
            names.append(rba_data.nmp_key(n))
        for n in rba_data.d_nucleotides:
            keys.append(rba_data.dntp_key(n))
            names.append(rba_data.dntp_key(n))        
        # key metabolites
        for m, name in rba_data.key_metabolites.items():
            keys.append(m)
            names.append(name)
        # cofactors
        for c in cofactors:
            keys.append(c.id)
            names.append(c.name)
        
        # retrieve items
        for key, name in zip(keys, names):
            # curated data available: nothing to do
            if key in self.metabolites: continue

            # no curated data: find metabolite id
            curated_data_added = True
            sbml_id = self._find(met_prefix + key + '_c', known_species)
            # use default concentration
            try:
                conc = rba_data.default_concentration[key]
            except KeyError:
                conc = 0
            self.metabolites[key] = Metabolite(name, sbml_id, conc)

        # write curated data to file (if necessary)
        if curated_data_added: self._write_curated_data()

    def _find(self, to_find, known_species):
        for s in known_species:
            if to_find.lower() == s.lower(): return s
        # not found:
        self._missing_information = True
        return self._missing_string

