from collections import defaultdict
import ConfigParser

import matplotlib.pyplot as plt
plt.style.use('ggplot')

from utils import *


def export_solution(ne_compartment_data, solution, config):
    fits_out_file = config.get('Output', 'fits_file')
    comp_out_file = config.get('Output', 'compartment_file')
    out_file_type = config.get('Output', 'file_type')

    if out_file_type == 'excel':
        solution.to_excel(fits_out_file)
        ne_compartment_data.to_excel(comp_out_file)
    elif out_file_type == 'csv':
        solution.to_csv(fits_out_file, sep=config.get('Output', 'delimiter'))
        ne_compartment_data.to_csv(comp_out_file, sep=config.get('Output', 'delimiter'))
    else:
        raise Exception('Output file type (%s) not allows. Use csv or excel.' % out_file_type)

def visualize_solution(exp_data, ratio_data, sol):

    growth_rates = np.array(exp_data['Growth rate'])
    comp2portion = defaultdict(list)
    for comp in ratio_data.index:
        for exp in exp_data.index:
            comp_p = ratio_data[exp][comp] / sum(ratio_data[exp])
            comp2portion[comp].append(comp_p)
    compartments = list(ratio_data.index)
    num_comp = len(compartments)
    col_num = 2
    row_num = int(np.ceil(num_comp / 2.))
    for i in range(num_comp):
        comp = compartments[i]
        exp_portions = ratio_data.loc[comp]
        intercept = sol.loc[comp]['intercept']
        slope = sol.loc[comp]['slope']
        fit = intercept + growth_rates * slope

        ax = plt.subplot(row_num, col_num, i+1)
        ax.plot(growth_rates, exp_portions, 'x')
        ax.plot(growth_rates, fit, '-')
        if i+1 in (4, 5):
            ax.set_xlabel('Growth rate (1/h)')
        ax.set_title(comp)
    plt.suptitle('Fraction of nonenzymatic protein')
    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    config = ConfigParser.ConfigParser()
    config.read('nonenz_per_compartment.cfg')

    experiment_data = load_experiment_data(config)
    protein_data = load_protein_data(config)
    location_data = load_location_data(config)
    category_data = load_category_data(config)
    comp_map_data = load_compartment_map(config)
    if comp_map_data is None:
        comp_map_data = create_mock_compartment_map(config, protein_data)
    is_enzymatic_data = load_is_enzymatic_data(config)

    # WE NEED ONLY THE PROTEINS THAT HAVE BOTH
    # LOCATION AND CATEGORY INFORMATION THAT
    # ARE NONENZYMATIC
    all_proteins = set(location_data.index) & set(category_data.index)
    cat_column = config.get('CategoryData', 'cat_column')
    is_enzymatic_col = config.get('IsEnzymaticData', 'is_enzymatic_column')
    
    access_array = []
    for prot in protein_data.index:
        if prot not in all_proteins:
            access_array.append(False)
            continue
        _cat = category_data.loc[prot][cat_column]
        if type(_cat) != unicode:
            _cat = _cat[0]
        if not is_enzymatic_data.loc[_cat][is_enzymatic_col]:
            access_array.append(True)
        else:
            access_array.append(False)

    ne_protein_data = protein_data[access_array]

    total_comp_data = compute_compartment_data(experiment_data, protein_data, 
                                               location_data, comp_map_data, config)
    ne_comp_data = compute_compartment_data(experiment_data, ne_protein_data, 
                                            location_data, comp_map_data, config)

    ratio_data = pd.DataFrame(index=total_comp_data.index, columns=total_comp_data.columns)
    for loc in total_comp_data.index:
        for exp in total_comp_data.columns:
            ratio_data.loc[loc, exp] = ne_comp_data.loc[loc, exp] / total_comp_data.loc[loc, exp]

    linear_fits = compute_linear_fits(experiment_data, ratio_data)

    export_solution(ne_comp_data, linear_fits, config)

    visualize_solution(experiment_data, ratio_data, linear_fits)

