import ConfigParser
import numpy as np

import matplotlib.pyplot as plt
plt.style.use('ggplot')

from utils import *

def get_exp_gr(exp_id, df, config):
    column_name = config.get('Experiment', 'growth_rate_column')
    return df.loc[exp_id][column_name]

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
    loc_data = load_location_data(config)
    if cm_data is None:
        cm_data = create_mock_compartment_map(config, prot_data)

    compartment_data = compute_compartment_data(exp_df, prot_data, loc_data, cm_data, config)

    solution = compute_linear_fits(exp_df, compartment_data, normalize=True, sums_to_one=True)

    export_solution(compartment_data, solution, config)

    visualize_solution(exp_df, compartment_data, solution)
    