from collections import defaultdict, MutableMapping
import numpy as np
import pandas as pd
from numpy.linalg import lstsq
from scipy.stats import linregress, pearsonr
from scipy.stats.mstats import gmean
import libsbml

import cplex
from cvxopt import solvers, matrix

R = 6.022140857e23

def numpy_to_cvxopt_matrix(A):
    if isinstance(A, np.ndarray):
        if A.ndim == 1:
            return matrix(A, (A.shape[0], 1), 'd')
        else:
            return matrix(A, A.shape, 'd')
    else:
        return A


def read_data(config, category):
    _file = config.get(category, 'file')
    _file_type = config.get(category, 'file_type')

    if _file_type == 'excel':
        sheet_name = config.get(category, 'sheet_name')
        skiprows = config.get(category, 'skiprows')
        df = pd.read_excel(_file,
                         sheet_name=sheet_name,
                         index_col=0,
                         skiprows=int(skiprows))
    elif _file_type == 'csv':
        delimiter = config.get(category, 'delimiter')
        df = pd.read_csv(_file,
                         index_col=0,
                         delimiter=delimiter)
    else:
        raise Exception('%s file type (%s) not allows. Use csv or excel.' % (category, exp_file_type))
    return df

def load_experiment_data(config):
    return read_data(config, 'Experiment')

def load_location_data(config):
    return read_data(config, 'LocationData')

def load_protein_data(config):
    return read_data(config, 'ProteinData')

def load_category_data(config):
    return read_data(config, 'CategoryData')

def load_compartment_map(config):
    if not config.has_section('CompartmentMap'):
        return None
    return read_data(config, 'CompartmentMap')

def load_is_enzymatic_data(config):
    return read_data(config, 'IsEnzymaticData')

def load_flux_data(config):
    return read_data(config, 'FluxData')


def create_mock_compartment_map(config, prot_data):
    cm_column = 'Mapped compartment'
    if not config.has_section('CompartmentMap'):
        config.add_section('CompartmentMap')
        config.set('CompartmentMap', 'cm_column', cm_column)
    location_column = config.get('ProteinData', 'location_column')
    locations = set([l.lower() for l in prot_data[location_column]])
    cm_data = pd.DataFrame(columns=[cm_column], index=locations)
    for loc in locations:
        cm_data.loc[loc][cm_column] = loc
    return cm_data

def compute_compartment_data(exp_df, prot_df, loc_df, cm_df, config):
    loc_col = config.get('LocationData', 'location_column')
    loc_df[loc_col] = loc_df[loc_col].str.lower()
    if loc_col not in prot_df.columns:
        total = pd.concat([prot_df, loc_df], join='inner', axis=1)
    else:
        prot_df[loc_col] = prot_df[loc_col].str.lower()
        total = prot_df
    locs = set(total[loc_col])

    cm_col = config.get('CompartmentMap', 'cm_column')
    
    import collections
    comp2comp = defaultdict(list)
    for c in cm_df.index:
        mapped = cm_df.loc[c][cm_col]
        comp2comp[mapped].append(c)

    comp_data = pd.DataFrame(columns=exp_df.index, index=comp2comp.keys())
    for loc in locs:
        try:
            true_loc = cm_df.loc[loc][cm_col]
        except Exception, e:
            import pdb
            pdb.set_trace()
        for exp in comp_data.columns:
            prot_w = sum(total[total[loc_col] == loc][exp])
            if np.isnan(comp_data[exp][true_loc]):
                comp_data[exp][true_loc] = prot_w
            else:
                comp_data[exp][true_loc] = comp_data[exp][true_loc] + prot_w
    return comp_data


def compute_linear_fits(exp_df, data_df, normalize=False, sums_to_one=False):
    sum_per_exp = []
    for exp in data_df.columns:
        sum_per_exp.append(sum(data_df[exp]))
    sum_per_exp = np.array(sum_per_exp)

    growth_rates = exp_df['Growth rate']
    comp2linfit = {}
    for comp in data_df.index:
        if normalize:
            vals = data_df.loc[comp] / sum_per_exp
        else:
            vals = data_df.loc[comp]
        slope, intercept, rval, pval, stderr = linregress(list(growth_rates), list(vals))
        if pval < 0.05:
            comp2linfit[comp] = True
        else:
            comp2linfit[comp] = False

    num_linfits = comp2linfit.values().count(True)
    num_constant = len(comp2linfit) - num_linfits
    # we don't fix the 'constant' compartments
    num_cols2 = 2 * len(comp2linfit)
    num_cols = 2 * num_linfits + num_constant
    num_exp = len(exp_df.index)
    if sums_to_one:
        num_rows = len(comp2linfit) * num_exp + num_exp
    else:
        num_rows = len(comp2linfit) * num_exp

    A = np.zeros((num_rows, num_cols))
    b = np.zeros((num_rows, 1))
    b_zero = np.zeros((num_rows, 1))
    gr = np.array(growth_rates)
    ones = np.ones((len(growth_rates),))

    idx_col = 0
    idx_row = 0
    last_constraint_idx = len(comp2linfit) * len(exp_df.index)

    for comp in data_df.index:
        if comp2linfit[comp]:
            A[idx_row: idx_row+num_exp, idx_col] = gr
            A[idx_row: idx_row+num_exp, idx_col + 1] = ones
            if sums_to_one:
                A[last_constraint_idx:, idx_col] = gr
                A[last_constraint_idx:, idx_col + 1] = ones
            if normalize:
                b[idx_row: idx_row+num_exp, 0] = np.array(data_df.loc[comp] / sum_per_exp)
            else:
                b[idx_row: idx_row+num_exp, 0] = np.array(data_df.loc[comp])
            idx_col = idx_col + 2
        else:
            if normalize:
                b[idx_row: idx_row+num_exp, 0] = np.array(data_df.loc[comp] / sum_per_exp)
            else:
                b[idx_row: idx_row+num_exp, 0] = np.array(data_df.loc[comp])
            A[idx_row: idx_row + num_exp, idx_col] = ones
            if sums_to_one:
                A[last_constraint_idx:, idx_col] = ones
            idx_col = idx_col + 1
        idx_row = idx_row + num_exp
    if sums_to_one:
        b[last_constraint_idx:, 0] = ones

    A2 = numpy_to_cvxopt_matrix(A)
    b2 = numpy_to_cvxopt_matrix(b)
    b_zero2 = numpy_to_cvxopt_matrix(b_zero)
    sol = solvers.qp(A2.T * A2, -(b2.T * A2).T, -A2, b_zero2)

    x, residuals, rank, s = lstsq(A, b)
    r, p = pearsonr(b, np.dot(A, x))
    r2, p2 = pearsonr(b2, np.dot(A2, sol['x']))
    print 'R^2 = ', r2[0]
    print 'p-value = ', p2[0]

    solution = pd.DataFrame(index=data_df.index, columns=['slope', 'intercept'])
    idx = 0
    x = sol['x'].T
    for comp in data_df.index:
        if comp2linfit[comp]:
            solution['slope'][comp] = x[idx]
            solution['intercept'][comp] = x[idx+1]
            idx = idx + 2
        else:
            solution['slope'][comp] = 0
            solution['intercept'][comp] = x[idx]
            idx = idx + 1
    return solution

def get_mmol_gcdw(protein_count, cell_dry_weight):
    return protein_count * 1e3 / (R * cell_dry_weight)

def get_gene_cnt_per_reaction(rba_model, fluxes=None):
    gene_cnt = defaultdict(int)
    for enz in rba_model.enzymes.enzymes:
        reaction_id = enz.id[2:-7]
        if 'duplicate' in reaction_id:
            dup_index = reaction_id.index('duplicate') - 1
            reaction_id = reaction_id[:dup_index]

        if fluxes is not None:
            if reaction_id in fluxes:
                for sr in enz.machinery_composition.reactants:
                    gene_cnt[sr.species] = gene_cnt[sr.species] + 1
        else:
            for sr in enz.machinery_composition.reactants:
                gene_cnt[sr.species] = gene_cnt[sr.species] + 1
    return gene_cnt

def is_enz_cytosolic(enz_id, rba_model):
    protein_compartments = defaultdict(int)
    enz = rba_model.enzymes.enzymes.get_by_id('R_%s_enzyme' % enz_id)
    for sp in enz.machinery_composition.reactants:
        prot = rba_model.proteins.macromolecules.get_by_id(sp.species)
        protein_compartments[prot.compartment] = protein_compartments[prot.compartment] + 1
    if len(protein_compartments) == 1 and protein_compartments.keys()[0] == 'Cytoplasm':
        return True
    else:
        return False

def create_fba_problem(sbml_file, biomass_reaction):
    # load model
    document = libsbml.readSBML(sbml_file)
    sbml_model = document.getModel()
    if sbml_model is None:
        raise Exception('Unable to load model from provided file path.')
    # load species and reactions
    reactions = {r[1].getId(): r[0] for r in enumerate(sbml_model.getListOfReactions())}
    if biomass_reaction not in reactions:
        raise Exception('Provided biomass reaction {} not in model {}.'.format(
            biomass_reaction,
            sbml_model.getName()
        ))
    species = {s[1].getId(): s[0] for s in enumerate(sbml_model.getListOfSpecies())}
    reaction_ids = sorted(reactions, key=reactions.get, reverse=False)
    species_ids = sorted(species, key=species.get, reverse=False)
    N_reac = len(reactions)
    N_spec = len(species)
    # create necessary matrices and vectors:
    # - stoichiometric matrix S
    # - upper / lower bound vectors: UB / LB
    # - right hand side vector: RHS
    # - objective function vector
    UB = np.zeros((N_reac,))
    LB = np.zeros((N_reac,))
    RHS = np.zeros((N_spec,))
    obj_fun_vec = np.zeros((N_reac,))
    obj_fun_vec[reactions[biomass_reaction]] = 1
    rows = defaultdict(list)
    values = defaultdict(list)

    for r_idx, reaction in enumerate(sbml_model.getListOfReactions()):
        plugin = reaction.getPlugin('fbc')
        if plugin is None:
            raise Exception('No fbc plugin defined for reaction {}. No way \
                to determine reaction bounds.'.format(reaction.getName()))
        reactants = {r.getSpecies(): r.getStoichiometry() for r in reaction.getListOfReactants()}
        products = {p.getSpecies(): p.getStoichiometry() for p in reaction.getListOfProducts()}
        reaction_species = list(set(reactants.keys() + products.keys()))
        stoichiometries = []
        for s in reaction_species:
            stoich = products.get(s, 0) - reactants.get(s, 0)
            stoichiometries.append(stoich)
        for sp, stoich in zip(reaction_species,stoichiometries):
            sp_idx = species[sp]
            rows[sp_idx].append(r_idx)
            values[sp_idx].append(stoich)

        UB[r_idx] = sbml_model.getParameter(plugin.getUpperFluxBound()).getValue()
        LB[r_idx] = sbml_model.getParameter(plugin.getLowerFluxBound()).getValue()
    # create CPLEX objects
    lp_problem = cplex.Cplex()
    cplex_rows = []
    for sp_idx in sorted(species.values()):
        cplex_rows.append(cplex.SparsePair(rows[sp_idx], values[sp_idx]))

    lp_problem.objective.set_sense(lp_problem.objective.sense.maximize)
    lp_problem.variables.add(names=reaction_ids, obj=obj_fun_vec,
                             ub=UB, lb=LB)
    lp_problem.linear_constraints.add(names=species_ids, rhs=RHS,
                                      senses=['E'] * len(species),
                                      lin_expr=cplex_rows)
    lp_problem.set_results_stream(None)
    return lp_problem

class Enzymes(MutableMapping):
    """A dictionary that applies an arbitrary key-altering
       function before accessing the keys"""

    def __init__(self, *args, **kwargs):
        self.store = dict()
        self.update(dict(*args, **kwargs))  # use the free update to set keys

    def __init__(self, rba_model, fluxes, protein_counts, protein_concentrations, gene_cnt, nz_gene_cnt, cdw):
        self.store = dict()
        self.fluxes = dict(fluxes)
        self.cdw = cdw
        for enz in rba_model.enzymes.enzymes:
            enz_id = enz.id[2:-7]
            is_cytosolic = is_enz_cytosolic(enz_id, rba_model)
            if 'duplicate' in enz_id:
                dup_index = enz_id.index('duplicate') - 1
                flux_id = enz_id[:dup_index]
            else:
                flux_id = enz_id
            #if flux_id not in fluxes:
            #    continue
            flux = 0. if flux_id not in fluxes else fluxes[flux_id]
            genes = {}
            for sr in enz.machinery_composition.reactants:
                gid = sr.species
                stoich = sr.stoichiometry
                r_cnt = gene_cnt[gid]
                nz_r_cnt = 0 if gid not in nz_gene_cnt else nz_gene_cnt[gid]
                count = 0 if gid not in protein_counts else protein_counts[gid]
                conc = 0 if gid not in protein_counts else protein_concentrations[gid]
                gene = Protein(gid, r_cnt, nz_r_cnt, stoich, count, conc)
                genes[gid] = gene
            enzyme = Enzyme(enz_id, flux, genes, is_cytosolic, cdw)

            self.store[enz_id] = enzyme
        self.correct_ATPS4rpp()
        self.adjust_isoenzymes()
        for enz in self.store.values():
            enz.set_kapp()

    def adjust_isoenzymes(self):
        flux2enzyme = defaultdict(list)
        for enz_id in self.store:
            if 'duplicate' in enz_id:
                dup_index = enz_id.index('duplicate') - 1
                store_id = enz_id[:dup_index]
            else:
                store_id = enz_id
            flux2enzyme[store_id].append(enz_id)
        for flux, enz_ids in flux2enzyme.items():
            enzymes = [self.store[_id] for _id in enz_ids]
            if len(enzymes) == 1:
                continue
            counts = np.array([enz.count for enz in enzymes])
            proportions = counts / np.sum(counts)
            for (enz, prop) in zip(enzymes, proportions):
                enz.flux = enz.flux * prop
        self.flux2enzyme = flux2enzyme

    def correct_ATPS4rpp(self):
        if 'ATPS4rpp' in self.store:
            ATPS4rpp_enz = self.store['ATPS4rpp']
            counts = np.array([gene.count for gene in ATPS4rpp_enz.genes.values()])
            stoichs = np.array([gene.stoichiometry for gene in ATPS4rpp_enz.genes.values()])
            reactions = np.array([gene.nz_reactions for gene in ATPS4rpp_enz.genes.values()])
            counts = counts / stoichs
            counts = counts / reactions
            counts = filter(lambda x: x != 0, counts)
            ATPS4rpp_enz.count = gmean(counts)


    def __getitem__(self, key):
        return self.store[self.__keytransform__(key)]

    def __setitem__(self, key, value):
        self.store[self.__keytransform__(key)] = value

    def __delitem__(self, key):
        del self.store[self.__keytransform__(key)]

    def __iter__(self):
        return iter(self.store)

    def __len__(self):
        return len(self.store)

    def __keytransform__(self, key):
        return key


class Enzyme(object):
    def __init__(self, name, flux, genes, is_cytosolic, cdw):
        self.name = name
        self.flux = flux
        self.genes = genes
        self.is_cytosolic = is_cytosolic
        self.cdw = cdw
        self.set_count()
        self.set_conc()
        self.set_kapp()


    def __repr__(self):
        _str = 'Reaction: {}\nFlux:    {}\nCount: {}\nkapp:      {}\n'.format(self.name, self.flux, self.count, self.kapp/3600)
        _str = _str + 'Genes\tReactions\tStoichiometry\tCount\n'
        gene_str = '\n'.join([str(gene) for gene in self.genes.values()])
        return _str + gene_str

    def set_count(self):
        counts = np.array([gene.count for gene in self.genes.values()])
        stoichs = np.array([gene.stoichiometry for gene in self.genes.values()])
        reactions = np.array([gene.nz_reactions for gene in self.genes.values()])
        # if there is any non-zero protein that is used only in this reaction, filter other zeros
        exclusive_non_zero = False
        if len(counts) > 1 and 0. in counts:
            for (cnt, reac) in zip(counts, reactions):
                if reac == 1 and cnt != 0:
                    exclusive_non_zero = True
        counts = counts / stoichs
        #counts = counts / reactions
        if exclusive_non_zero:
            counts = filter(lambda x: x != 0, counts)
        self.count = gmean(counts)

    def set_conc(self):
        self.set_count()
        self.concentration = self.count / R / self.cdw * 1e3

    def set_kapp(self):
        self.kapp = self.flux / self.concentration if (self.concentration != 0 and not np.isnan(self.concentration)) else 0


class Protein(object):
    def __init__(self, gene, reactions, nz_reactions, stoichiometry, count, concentration):
        self.gene = gene
        self.reactions = reactions
        self.nz_reactions = nz_reactions
        self.stoichiometry = stoichiometry
        self.count = count
        self.concentration = concentration

    def __repr__(self):
        _str = '{}\t{}/{}\t{}\t{}\t{}'.format(self.gene,
                                           self.nz_reactions,
                                           self.reactions,
                                           self.stoichiometry,
                                           self.count,
                                           self.concentration)
        return _str

    def __str__(self):
        return self.__repr__()
