from collections import defaultdict
import numpy as np
import pandas as pd
from numpy.linalg import lstsq
from scipy.stats import linregress, pearsonr

from cvxopt import solvers, matrix


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
		if comp == 'secreted':
			comp2linfit[comp] = False
		elif pval < 0.05:
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