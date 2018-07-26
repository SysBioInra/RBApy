from collections import defaultdict
import ConfigParser
import numpy as np
import pandas as pd
from numpy.linalg import lstsq
from scipy.stats import pearsonr
from scipy.stats import linregress
import matplotlib.pyplot as plt
plt.style.use('ggplot')

from cvxopt import solvers, matrix


def numpy_to_cvxopt_matrix(A):
    if isinstance(A, np.ndarray):
        if A.ndim == 1:
            return matrix(A, (A.shape[0], 1), 'd')
        else:
            return matrix(A, A.shape, 'd')
    else:
        return A



def get_exp_gr(exp_id, df, config):
	column_name = config.get('Experiment', 'growth_rate_column')
	return df.loc[exp_id][column_name]

def load_experiment_data(config):
	exp_file = config.get('Experiment', 'file')
	exp_file_type = config.get('Experiment', 'file_type')

	if exp_file_type == 'excel':
		sheet_name = config.get('Experiment', 'sheet_name')
		skiprows = config.get('Experiment', 'skiprows')
		df = pd.read_excel(exp_file,
						 sheet_name=sheet_name,
						 index_col=0,
						 skiprows=int(skiprows))
	elif exp_file_type == 'csv':
		delimiter = config.get('Experiment', 'delimiter')
		df = pd.read_csv(exp_file,
						 index_col=0,
						 delimiter=delimiter)
	else:
		raise Exception('Experiment file type (%s) not allows. Use csv or excel.' % exp_file_type)
	return df

def load_protein_data(config):
	prot_file = config.get('ProteinData', 'file')
	prot_file_type = config.get('ProteinData', 'file_type')
	if prot_file_type == 'excel':
		sheet_name = config.get('ProteinData', 'sheet_name')
		skiprows = config.get('ProteinData', 'skiprows')
		df = pd.read_excel(prot_file,
						 sheet_name=sheet_name,
						 index_col=0,
						 skiprows=int(skiprows))
	elif prot_file_type == 'csv':
		delimiter = config.get('ProteinData', 'delimiter')
		df = pd.read_csv(prot_file,
						 index_col=0,
						 delimiter=delimiter)
	else:
		raise Exception('ProteinData file type (%s) not allows. Use csv or excel.' % prot_file_type)
	return df

def load_compartment_map(config):
	try:
		cm_file = config.get('CompartmentMap', 'file')
	except Exception, e:
		return None
	cm_file_type = config.get('CompartmentMap', 'file_type')
	if cm_file_type == 'excel':
		sheet_name = config.get('CompartmentMap', 'sheet_name')
		skiprows = config.get('CompartmentMap', 'skiprows')
		df = pd.read_excel(cm_file,
						 sheet_name=sheet_name,
						 index_col=0,
						 skiprows=int(skiprows))
	elif cm_file_type == 'csv':
		delimiter = config.get('CompartmentMap', 'delimiter')
		df = pd.read_csv(cm_file,
						 index_col=0,
						 delimiter=delimiter)
	else:
		raise Exception('CompartmentMap file type (%s) not allows. Use csv or excel.' % cm_file_type)
	return df

def compute_compartment_data(exp_df, prot_df, cm_df, config):
	loc_col = config.get('ProteinData', 'location_column')
	cm_col = config.get('CompartmentMap', 'cm_column')
	prot_df[loc_col] = prot_df[loc_col].str.lower()
	locs = set(prot_df[loc_col])
	
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
			prot_w = sum(prot_df[prot_df[loc_col] == loc][exp])
			if np.isnan(comp_data[exp][true_loc]):
				comp_data[exp][true_loc] = prot_w
			else:
				comp_data[exp][true_loc] = comp_data[exp][true_loc] + prot_w
	return comp_data

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



def compute_linear_fits(exp_df, compartment_data):
	sum_per_exp = []
	for exp in compartment_data.columns:
		sum_per_exp.append(sum(compartment_data[exp]))
	sum_per_exp = np.array(sum_per_exp)

	growth_rates = exp_df['Growth rate']
	comp2linfit = {}
	for comp in compartment_data.index:
		vals = compartment_data.loc[comp] / sum_per_exp
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
	num_rows = len(comp2linfit) * num_exp + num_exp

	A = np.zeros((num_rows, num_cols))
	b = np.zeros((num_rows, 1))
	b_zero = np.zeros((num_rows, 1))
	gr = np.array(growth_rates)
	ones = np.ones((len(growth_rates),))

	idx_col = 0
	idx_row = 0
	last_constraint_idx = len(comp2linfit) * len(exp_df.index)

	for comp in compartment_data.index:
		if comp2linfit[comp]:
			A[idx_row: idx_row+num_exp, idx_col] = gr
			A[idx_row: idx_row+num_exp, idx_col + 1] = ones
			A[last_constraint_idx:, idx_col] = gr
			A[last_constraint_idx:, idx_col + 1] = ones
			b[idx_row: idx_row+num_exp, 0] = np.array(compartment_data.loc[comp] / sum_per_exp)
			idx_col = idx_col + 2
		else:
			b[idx_row: idx_row+num_exp, 0] = np.array(compartment_data.loc[comp] / sum_per_exp)
			A[idx_row: idx_row + num_exp, idx_col] = ones
			A[last_constraint_idx:, idx_col] = ones
			idx_col = idx_col + 1
		idx_row = idx_row + num_exp
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

	solution = pd.DataFrame(index=compartment_data.index, columns=['slope', 'intercept'])
	idx = 0
	x = sol['x'].T
	for comp in compartment_data.index:
		if comp2linfit[comp]:
			solution['slope'][comp] = x[idx]
			solution['intercept'][comp] = x[idx+1]
			idx = idx + 2
		else:
			solution['slope'][comp] = 0
			solution['intercept'][comp] = x[idx]
			idx = idx + 1
	return solution

def export_solution(compartment_data, solution, config):
	fits_out_file = config.get('Output', 'fits_file')
	comp_out_file = config.get('Output', 'compartment_file')
	out_file_type = config.get('Output', 'file_type')

	if out_file_type == 'excel':
		solution.to_excel(fits_out_file)
		compartment_data.to_excel(comp_out_file)
	elif out_file_type == 'csv':
		solution.to_csv(fits_out_file, sep=config.get('Output', 'delimiter'))
		compartment_data.to_csv(comp_out_file, sep=config.get('Output', 'delimiter'))
	else:
		raise Exception('Output file type (%s) not allows. Use csv or excel.' % out_file_type)

def visualize_solution(exp_data, comp_data, sol):

	growth_rates = np.array(exp_data['Growth rate'])
	comp2portion = defaultdict(list)
	for comp in comp_data.index:
		for exp in exp_data.index:
			comp_p = comp_data[exp][comp] / sum(comp_data[exp])
			comp2portion[comp].append(comp_p)


	compartments = list(comp_data.index)
	num_comp = len(compartments)
	col_num = 2
	row_num = int(np.ceil(num_comp / 2.))
	for i in range(num_comp):

		comp = compartments[i]
		exp_portions = comp2portion[comp]
		intercept = sol.loc[comp]['intercept']
		slope = sol.loc[comp]['slope']
		fit = intercept + growth_rates * slope

		ax = plt.subplot(row_num, col_num, i+1)
		ax.plot(growth_rates, exp_portions, 'x')
		ax.plot(growth_rates, fit, '-')
		if i+1 in (4, 5):
			ax.set_xlabel('Growth rate (1/h)')
		ax.set_title(comp)

	plt.suptitle('Fraction of protein per compartment')
	plt.tight_layout()
	plt.show()


if __name__ == '__main__':
	config = ConfigParser.ConfigParser()
	config.read('prot_per_compartment.cfg')

	exp_df = load_experiment_data(config)
	prot_data = load_protein_data(config)
	cm_data = load_compartment_map(config)
	if cm_data is None:
		cm_data = create_mock_compartment_map(config, prot_data)
	compartment_data = compute_compartment_data(exp_df, prot_data, cm_data, config)
	solution = compute_linear_fits(exp_df, compartment_data)
	export_solution(compartment_data, solution, config)
	visualize_solution(exp_df, compartment_data, solution)
	