import ConfigParser
import pandas as pd

import rba

from utils import *

import cplex


def set_bounds(fba_model, flux_data, config):
    lb_column = config.get('FluxData', 'lb_column')
    ub_column = config.get('FluxData', 'ub_column')
    variable_names = fba_model.variables.get_names()
    for fluxes in flux_data.index:
        reactions = fluxes.split(',')
        indices = [variable_names.index(reac) for reac in reactions]
        lin_expr = [cplex.SparsePair(indices, [1]*len(indices)),]

        names_lb = ['{}_LB'.format(fluxes),]
        senses_lb = ['GE',]
        rhs_lb = [flux_data[lb_column][fluxes],]
        fba_model.linear_constraints.add(names=names_lb)
        fba_model.linear_constraints.set_linear_components(zip(names_lb, lin_expr))
        fba_model.linear_constraints.set_senses(zip(names_lb, senses_lb))
        fba_model.linear_constraints.set_rhs(zip(names_lb, rhs_lb))

        names_ub = ['{}_UB'.format(fluxes),]
        senses_ub = ['LE',]
        rhs_ub = [flux_data[ub_column][fluxes],]
        fba_model.linear_constraints.add(names=names_ub)
        fba_model.linear_constraints.set_linear_components(zip(names_ub, lin_expr))
        fba_model.linear_constraints.set_senses(zip(names_ub, senses_ub))
        fba_model.linear_constraints.set_rhs(zip(names_ub, rhs_ub))
    return fba_model


if __name__ == '__main__':
    config = ConfigParser.ConfigParser()
    config.read('kapp.cfg')

    flux_data = load_flux_data(config)

    protein_data = load_protein_data(config)
    count_column = config.get('ProteinData', 'count_column')
    cdw = float(config.get('ProteinData', 'cdw'))

    fba_model_file = config.get('Model', 'fba')
    biomass_reaction = config.get('Model', 'biomass_reaction')
    fba_model = create_fba_problem(fba_model_file, biomass_reaction)
    fba_model = set_bounds(fba_model, flux_data, config)
    fba_model.solve()
    fluxes = {a[2:]: b for a,b in zip(fba_model.variables.get_names(),
                                  fba_model.solution.get_values())}
    non_zero_fluxes = dict(filter(lambda i : i[1] != 0, fluxes.items()))

    rba_model_dir = config.get('Model', 'rba')
    rba_model = rba.RbaModel.from_xml(rba_model_dir)

    gene_cnt = get_gene_cnt_per_reaction(rba_model)
    nz_gene_cnt = get_gene_cnt_per_reaction(rba_model, non_zero_fluxes)

    counts = {p: protein_data.loc[p][count_column] for p in protein_data.index}
    concs = {p[0]: get_mmol_gcdw(p[1], cdw) for p in counts.items()}

    enzymes = Enzymes(rba_model, non_zero_fluxes, counts, concs, gene_cnt, nz_gene_cnt, cdw)

    enz_ids = [enz for enz in enzymes.keys() if enzymes[enz].is_cytosolic]
    enz_kapp_df = pd.DataFrame(index=enz_ids, columns=['kapp'])
    enz_kapp_df.index.name = 'Reaction ID'

    for enz_id in enz_ids:
        enz = enzymes[enz_id]
        enz_kapp_df.loc[enz_id]['kapp'] = enz.kapp / 3600

    output_file = config.get('Output', 'file')
    output_file_type = config.get('Output', 'file_type')
    nz_values = enz_kapp_df.iloc[enz_kapp_df['kapp'].nonzero()]
    if output_file_type == 'excel':
        nz_values.to_excel(output_file)
    elif output_file_type == 'csv':
        delimiter = config.get('Output', 'delimiter')
        nz_values.to_csv(output_file, str(delimiter))
