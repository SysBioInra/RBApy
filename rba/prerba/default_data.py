"""Module storing data and data-related functions."""

# python 2/3 compatibility
from __future__ import division, print_function, absolute_import

# local imports
import rba.xml
from rba.xml import Function

GROWTH_RATE = 'growth_rate'


def build_aggregate(id_, fn_refs):
    """Build aggregate with given identifiers and function references."""
    result = rba.xml.Aggregate(id_, 'multiplication')
    for ref in fn_refs:
        result.function_references.append(rba.xml.FunctionReference(ref))
    return result


class DefaultData(object):
    """Class holding default RBA data."""

    def __init__(self):
        self.parameters = DefaultParameters()
        self.metabolites = DefaultMetabolites()
        self.activity = DefaultActivity()
        self.atpm_reaction = 'R_maintenance_atp'


class DefaultParameters(object):
    """Class holding default RBA parameter data."""

    def metabolite_concentration(self, id_):
        return id_ + '_concentration'

    def metabolite_concentration_function(self, id_, concentration):
        return Function(self.metabolite_concentration(id_),
                        'constant', {'CONSTANT': concentration})

    def all_functions(self):
        pass

    def process_functions(self):
        return [
            # zero function
            Function('zero', 'constant', {'CONSTANT': 0}),
            # efficiencies
            Function('ribosome_efficiency_MM', 'michaelisMenten',
                     {'kmax': 97200, 'Km': 0.5, 'Y_MIN': 32400}),
            Function('ribosome_efficiency_CM', 'constant',
                     {'CONSTANT': 97200}),
            Function('fraction_active_ribosomes', 'exponential',
                     {'RATE': -0.083333}),
            Function('chaperone_efficiency_LM', 'linear',
                     {'LINEAR_COEF': 36044.48, 'LINEAR_CONSTANT': -2888.0051,
                      'X_MIN': 0.25, 'X_MAX': 1.6,
                      'Y_MIN': float('-Inf'), 'Y_MAX': float('Inf')}),
            # targets
            Function('mrna_degradation_flux', 'constant',
                     {'CONSTANT': 0.15996}),
            Function('mrna_concentration', 'constant', {'CONSTANT': 0.01}),
            Function('dna_concentration', 'constant', {'CONSTANT': 0.0807}),
            Function('maintenance_atp', 'linear',
                     {'LINEAR_COEF': 12.1595, 'LINEAR_CONSTANT': -3.1595,
                      'X_MIN': 1, 'X_MAX': float('Inf'),
                      'Y_MIN': float('-Inf'), 'Y_MAX': float('Inf')})
            ]

    def process_aggregates(self):
        return [build_aggregate(
            'ribosome_capacity',
            ['ribosome_efficiency_MM', 'fraction_active_ribosomes']
            )]

    def density_functions(self, cytoplasm, external, other):
        cytoplasm_density = 4.8972
        cytoplasm_fraction = 0.7
        external_fraction = 0.1
        other_fraction = 0.2 / len(other)
        fns = [Function('amino_acid_concentration', 'linear',
                        {'LINEAR_COEF': -0.9757, 'LINEAR_CONSTANT': 6.3138,
                         'X_MIN': 0.25, 'X_MAX': 1.6,
                         'Y_MIN': float('-Inf'), 'Y_MAX': float('Inf')})]
        # protein fraction
        fns.append(self.protein_fraction(cytoplasm, cytoplasm_fraction))
        fns.append(self.protein_fraction(external, external_fraction))
        fns += [self.protein_fraction(c, other_fraction) for c in other]
        # non enzymatic fraction
        fns.append(self.non_enzymatic_fraction_cytoplasm(cytoplasm))
        fns.append(self.non_enzymatic_fraction_secreted(external))
        fns += [self.non_enzymatic_fraction_other(c) for c in other]
        # total density
        fns.append(Function(cytoplasm + '_density',
                            'constant', {'CONSTANT': cytoplasm_density}))
        return fns

    def density_aggregates(self, cytoplasm, external, other):
        aggregates = []
        for cpt in other:
            aggregates.append(build_aggregate(
                cpt + '_density',
                ['amino_acid_concentration', self.protein_fraction_id(cpt)]
                ))
        # non enzymatic density
        for cpt in [cytoplasm, external] + other:
            aggregates.append(build_aggregate(
                'nonenzymatic_proteins_' + cpt,
                ['amino_acid_concentration',
                 'inverse_average_protein_length',
                 self.protein_fraction_id(cpt),
                 self.non_enzymatic_fraction_id(cpt)]
                ))
        return aggregates

    @staticmethod
    def inverse_average_protein_length(length):
        """
        Return xml structure containing average protein length.
        """
        return Function('inverse_average_protein_length', 'constant',
                        {'CONSTANT': 1.0/length})

    @staticmethod
    def protein_fraction_id(compartment_id):
        """
        Return function identifier for protein fraction in given compartment.
        """
        return 'fraction_protein_' + compartment_id

    def protein_fraction(self, compartment_id, fraction):
        """
        Return xml structure containing protein fraction in given compartment.
        """
        return Function(self.protein_fraction_id(compartment_id),
                        'constant', {'CONSTANT': fraction})

    @staticmethod
    def non_enzymatic_fraction_id(compartment_id):
        """
        Return function identifier for non-enzymatic fraction.
        """
        return 'fraction_non_enzymatic_protein_' + compartment_id

    def non_enzymatic_fraction_cytoplasm(self, compartment_id):
        """
        Return xml structure containing non-enzymatic fraction in cytosol.
        """
        params = {'LINEAR_COEF': -0.0502, 'LINEAR_CONSTANT': 0.3526,
                  'X_MIN': 0.25, 'X_MAX': 1.6, 'Y_MIN': 0, 'Y_MAX': 1}
        return Function(self.non_enzymatic_fraction_id(compartment_id),
                        'linear', params)

    def non_enzymatic_fraction_secreted(self, compartment_id):
        """
        Return xml structure containing external non-enzymatic fraction.
        """
        return Function(self.non_enzymatic_fraction_id(compartment_id),
                        'constant', {'CONSTANT': 1})

    def non_enzymatic_fraction_other(self, compartment_id):
        """
        Return xml structure containing generic non-enzymatic fraction.
        """
        params = {'LINEAR_COEF': -0.0496, 'LINEAR_CONSTANT': 0.8965,
                  'X_MIN': 0.25, 'X_MAX': 1.6, 'Y_MIN': 0, 'Y_MAX': 1}
        return Function(self.non_enzymatic_fraction_id(compartment_id),
                        'linear', params)


class DefaultMetabolites(object):
    """
    Class holding default metabolite data.
    """

    def __init__(self):
        # key metabolites and components
        self.aas_3L = {
            'A': 'ala', 'C': 'cys', 'D': 'asp', 'E': 'glu', 'F': 'phe',
            'G': 'gly', 'H': 'his', 'I': 'ile', 'K': 'lys', 'L': 'leu',
            'fM': 'fmet', 'M': 'met', 'N': 'asn', 'P': 'pro', 'Q': 'gln',
            'R': 'arg', 'S': 'ser', 'T': 'thr', 'V': 'val', 'W': 'trp',
            'Y': 'tyr'
            }
        self.aas = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M',
                    'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
        self.aa_fM = 'fM'
        self.nucleotides = ['A', 'C', 'G', 'U']
        self.d_nucleotides = ['A', 'C', 'G', 'T']
        self.key_metabolites = {
            'Pi': 'Phosphate', 'COA': 'Coenzyme A', 'NAD': 'NAD',
            'NADP': 'NADP', 'NADPH': 'NADPH', 'K': 'Potassium',
            'SPMD': 'SPMD', 'H2O': 'H2O', 'H': 'Proton',
            'PPi': 'Pyrophosphate',
            '10FTHF': '10-Formyltetrahydrofolate', 'FOR': 'Formate'
            }
        self.mrna = 'mrna'
        self.dna = 'dna'
        # default concentrations
        self.concentration = {
            'COA': 0.0003, 'NAD': 0.0162, 'NADP': 0.0009, 'NADPH': 0.0002,
            '10FTHF': 0.0004, 'K': 0.7063, 'SPMD': 0.007,
            'Pi': 0.0144, 'PPi': 0.0009,
            'AMP': 0.0047, 'CMP': 0.001, 'GMP': 0.0005, 'UMP': 0,
            'ADP': 0.0026, 'CDP': 0.0003, 'GDP': 0.0002, 'UDP': 0,
            'ATP': 0.003, 'CTP': 0.0005, 'GTP': 0.0004, 'UTP': 0
            }
        # trna concentration is actually a nucleotide concentration
        trna_concentration = {
            'C': 0.000642, 'P': 0.001550, 'H': 0.000664, 'D': 0.001690,
            'S': 0.001513, 'Q': 0.001845, 'I': 0.002037, 'M': 0.001078,
            'K': 0.002406, 'T': 0.001779, 'F': 0.001299, 'A': 0.003602,
            'G': 0.004295, 'E': 0.001845, 'L': 0.003159, 'R': 0.002074,
            'W': 0.000399, 'N': 0.001690, 'Y': 0.000967, 'V': 0.002967
            }
        av_trna_length = 75
        for trna, conc in trna_concentration.items():
            self.concentration[self.uncharged_trna_key(trna)] \
                = conc / av_trna_length

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

    @staticmethod
    def average_protein_id(compartment_id):
        """
        Return identifier of average protein in given compartment.
        """
        return 'average_protein_' + compartment_id

    def process_metabolites(self):
        """
        Return internal ids and names of metabolites involved in processes.
        """
        # metabolites to retrieve
        keys = []  # internal id of metabolite
        names = []  # name of metabolite
        # methionine
        keys.append('MET')
        names.append('Methionine')
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


class DefaultActivity(object):
    """Class holding default RBA enzyme activity data."""

    def __init__(self):
        """
        Default constructor.
        """
        # default ids
        self.efficiency_id = 'default_efficiency'
        self.transport_id = 'default_transporter_efficiency'
        # default medium
        self.medium_concentration = 10

    def efficiency_function(self):
        return Function(self.efficiency_id, 'constant',
                        {'CONSTANT': 200000}, GROWTH_RATE)

    def transport_function(self):
        return Function(self.transport_id, 'constant',
                        {'CONSTANT': 2e6}, GROWTH_RATE)

    def transport_aggregate_id(self, reaction):
        return '{}_efficiency'.format(reaction)

    def transport_functions(self, reaction, metabolites):
        return [
            Function(self._transport_fn_id(reaction, met),
                     'michaelisMenten', {'Km': 0.8, 'kmax': 1}, met)
            for met in metabolites
        ]

    def _transport_fn_id(self, reaction, metabolite):
        return '{}_{}_transport_factor'.format(reaction, metabolite)

    def transport_aggregate(self, reaction, metabolites):
        agg_fns = ([self.transport_id]
                   + [self._transport_fn_id(reaction, m) for m in metabolites])
        return build_aggregate(self.transport_aggregate_id(reaction), agg_fns)
