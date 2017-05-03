
from ..rba_xml import Function

## function information
cytoplasm_density = 4.8972
cytoplasm_fraction = 0.7
secreted_fraction = 0.1
other_fraction = 0.2

def inverse_average_protein_length(length):
    return Function('inverse_average_protein_length', 'constant',
                    {'CONSTANT': 1.0/length})

def protein_fraction_id(compartment_id):
    return 'fraction_protein_' + compartment_id

def protein_fraction(compartment_id, fraction):
    return Function(protein_fraction_id(compartment_id),
                    'constant', {'CONSTANT': fraction})

def non_enzymatic_fraction_id(compartment_id):
    return 'fraction_non_enzymatic_protein_' + compartment_id

def non_enzymatic_fraction_cytoplasm(compartment_id):
    return Function(non_enzymatic_fraction_id(compartment_id), 'linear',
                    {'LINEAR_COEF': -0.0502, 'LINEAR_CONSTANT': 0.3526,
                     'X_MIN': 0.25, 'X_MAX': 1.6, 'Y_MIN': 0, 'Y_MAX': 1})

def non_enzymatic_fraction_secreted(compartment_id):
    return Function(non_enzymatic_fraction_id(compartment_id), 'constant',
                    {'CONSTANT': 1})

def non_enzymatic_fraction_other(compartment_id):
    return Function(non_enzymatic_fraction_id(compartment_id), 'linear',
                    {'LINEAR_COEF': -0.0496, 'LINEAR_CONSTANT': 0.8965,
                     'X_MIN': 0.25, 'X_MAX': 1.6, 'Y_MIN': 0, 'Y_MAX': 1})

fns = []
fns.append(Function('amino_acid_concentration', 'linear',
                    {'LINEAR_COEF': -0.9757, 'LINEAR_CONSTANT': 6.3138,
                     'X_MIN': 0.25, 'X_MAX': 1.6,
                     'Y_MIN': float('-Inf'), 'Y_MAX': float('Inf')}))
fns.append(Function('flagella_speed', 'constant', {'CONSTANT': 5.81}))
fns.append(Function('flagella_h_consumption', 'constant', {'CONSTANT': 0.9415}))
fns.append(Function('number_flagella', 'linear',
                    {'LINEAR_COEF': 4.5197, 'LINEAR_CONSTANT': 3.7991,
                     'X_MIN': 0.25, 'X_MAX': 1.6,
                     'Y_MIN': float('-Inf'), 'Y_MAX': float('Inf')}))
fns.append(Function('ribosome_efficiency_MM', 'michaelisMenten',
                    {'kmax': 97200, 'Km': 0.5, 'Y_MIN': 32400}))
fns.append(Function('ribosome_efficiency_CM', 'constant', {'CONSTANT': 97200}))
fns.append(Function('fraction_active_ribosomes', 'exponential',
                    {'RATE': -0.083333}))
fns.append(Function('chaperone_efficiency_LM', 'linear',
                    {'LINEAR_COEF': 36044.48, 'LINEAR_CONSTANT': -2888.0051,
                     'X_MIN': 0.25, 'X_MAX': 1.6,
                     'Y_MIN': float('-Inf'), 'Y_MAX': float('Inf')}))
fns.append(Function('maintenance_atp', 'linear',
                    {'LINEAR_COEF': 12.1595, 'LINEAR_CONSTANT': -3.1595,
                     'X_MIN': 1, 'X_MAX': float('Inf'),
                     'Y_MIN': float('-Inf'), 'Y_MAX': float('Inf')}))

# additional reactions
atpm_reaction = 'R_maintenance_atp'

# key metabolites and components
aas_3L = {'A':'ala', 'C':'cys', 'D':'asp', 'E':'glu', 'F':'phe', 'G':'gly',
          'H':'his', 'I':'ile', 'K':'lys', 'L':'leu', 'fM': 'fmet','M':'met',
          'N':'asn', 'P':'pro', 'Q':'gln', 'R':'arg', 'S':'ser', 'T':'thr',
          'V':'val', 'W':'trp', 'Y':'tyr'}
aas = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R',
       'S','T','V','W','Y']
aa_fM = 'fM'
nucleotides = ['A', 'C', 'G', 'U']
d_nucleotides = ['A', 'C', 'G', 'T']
key_metabolites = {'Pi': 'Phosphate', 'COA': 'Coenzyme A', 'NAD': 'NAD',
                   'NADP': 'NADP', 'NADPH': 'NADPH', 'K': 'Potassium',
                   'SPMD': 'SPMD', 'H2O': 'H2O', 'H': 'Proton',
                   'PPi': 'Pyrophosphate',
                   '10FTHF': '10-Formyltetrahydrofolate',
                   'FOR': 'Formate'}
mrna = 'mrna'
dna = 'dna'

def charged_trna_key(aa):
    return aas_3L[aa].upper() + 'TRNA'

def charged_trna_name(aa):
    return 'Charged trna ' + aas_3L[aa]

def uncharged_trna_key(aa):
    return 'TRNA' + aas_3L[aa].upper()

def uncharged_trna_name(aa):
    return 'Uncharged trna ' + aas_3L[aa]

def ntp_key(n):
    return n + 'TP'

def ndp_key(n):
    return n + 'DP'

def nmp_key(n):
    return n + 'MP'

def dntp_key(n):
    return 'd' + n + 'TP'

def average_protein_id(compartment_id):
    return 'average_protein_' + compartment_id

## default concentrations
default_concentration \
    = {'COA': 0.0003, 'NAD': 0.0162, 'NADP':0.0009, 'NADPH': 0.0002,
       '10FTHF': 0.0004, 'K': 0.7063, 'SPMD': 0.007,
       'Pi': 0.0144, 'PPi': 0.0009,
       'AMP': 0.0047, 'CMP': 0.001, 'GMP': 0.0005, 'UMP': 0,
       'ADP': 0.0026, 'CDP': 0.0003, 'GDP': 0.0002, 'UDP': 0,
       'ATP': 0.003, 'CTP': 0.0005, 'GTP': 0.0004, 'UTP': 0}
# trna concentration is actually a nucleotide concentration
trna_concentration \
    = { 'C': 0.000642, 'P': 0.001550, 'H': 0.000664, 'D': 0.001690,
        'S': 0.001513, 'Q': 0.001845, 'I': 0.002037, 'M': 0.001078,
        'K': 0.002406, 'T': 0.001779, 'F': 0.001299, 'A': 0.003602,
        'G': 0.004295, 'E': 0.001845, 'L': 0.003159, 'R': 0.002074,
        'W': 0.000399, 'N': 0.001690, 'Y': 0.000967, 'V': 0.002967 }
av_trna_length = 75
for trna, conc in trna_concentration.iteritems():
    default_concentration[uncharged_trna_key(trna)] = conc / av_trna_length

## default catalytic activity
default_catalytic_activity = 200000
default_transporter_activity = 2e6
default_import_Km = 0.8
default_import_kmax = 1

## default medium
default_medium_concentration = 10
