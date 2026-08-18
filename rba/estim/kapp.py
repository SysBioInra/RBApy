import configparser
import pandas as pd

import rba

from utils import *

import cplex

def set_bounds_parsimonious_fba(fba_model, flux_data, config):
    lb_column = config.get('FluxData', 'lb_column')
    ub_column = config.get('FluxData', 'ub_column')
    variable_names = fba_model.variables.get_names()
    for fluxes in flux_data.index:
        reactions = fluxes.split(',')
        if float(flux_data[ub_column][fluxes]) >= 0:
            if "{}__FW".format(fluxes) in variable_names:
                parsimonious_fba_rxn_for_ub=str("{}__FW".format(fluxes))
            else:
                parsimonious_fba_rxn_for_ub=fluxes
        elif float(flux_data[ub_column][fluxes]) < 0:
            if "{}__BW".format(fluxes) in variable_names:
                parsimonious_fba_rxn_for_ub=str("{}__BW".format(fluxes))
            else:
                parsimonious_fba_rxn_for_ub=fluxes
        if float(flux_data[lb_column][fluxes]) >= 0:
            if "{}__FW".format(fluxes) in variable_names:
                parsimonious_fba_rxn_for_lb=str("{}__FW".format(fluxes))
            else:
                parsimonious_fba_rxn_for_ub=fluxes
        elif float(flux_data[lb_column][fluxes]) < 0:
            if "{}__BW".format(fluxes) in variable_names:
                parsimonious_fba_rxn_for_lb=str("{}__BW".format(fluxes))
            else:
                parsimonious_fba_rxn_for_ub=fluxes

        names_lb = ['{}_LB'.format(parsimonious_fba_rxn_for_lb),]
        senses_lb = ['GE',]
        rhs_lb = [float(flux_data[lb_column][fluxes]),]
        fba_model.linear_constraints.add(names=names_lb)
        fba_model.linear_constraints.set_linear_components(zip(names_lb, [cplex.SparsePair([variable_names.index(parsimonious_fba_rxn_for_lb)], [1.0]),]))
        fba_model.linear_constraints.set_senses(zip(names_lb, senses_lb))
        fba_model.linear_constraints.set_rhs(zip(names_lb, rhs_lb))

        names_ub = ['{}_UB'.format(parsimonious_fba_rxn_for_ub),]
        senses_ub = ['LE',]
        rhs_ub = [float(flux_data[ub_column][fluxes]),]
        fba_model.linear_constraints.add(names=names_ub)
        fba_model.linear_constraints.set_linear_components(zip(names_ub, [cplex.SparsePair([variable_names.index(parsimonious_fba_rxn_for_ub)], [1.0]),]))
        fba_model.linear_constraints.set_senses(zip(names_ub, senses_ub))
        fba_model.linear_constraints.set_rhs(zip(names_ub, rhs_ub))
    return fba_model

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
        rhs_lb = [float(flux_data[lb_column][fluxes]),]
        fba_model.linear_constraints.add(names=names_lb)
        fba_model.linear_constraints.set_linear_components(zip(names_lb, lin_expr))
        fba_model.linear_constraints.set_senses(zip(names_lb, senses_lb))
        fba_model.linear_constraints.set_rhs(zip(names_lb, rhs_lb))

        names_ub = ['{}_UB'.format(fluxes),]
        senses_ub = ['LE',]
        rhs_ub = [float(flux_data[ub_column][fluxes]),]
        fba_model.linear_constraints.add(names=names_ub)
        fba_model.linear_constraints.set_linear_components(zip(names_ub, lin_expr))
        fba_model.linear_constraints.set_senses(zip(names_ub, senses_ub))
        fba_model.linear_constraints.set_rhs(zip(names_ub, rhs_ub))
    return fba_model


if __name__ == '__main__':
    config = configparser.ConfigParser()
    config.read('kapp.cfg')

    flux_data = load_flux_data(config)

    protein_data = load_protein_data(config)
    count_column = config.get('ProteinData', 'count_column')
    cdw = float(config.get('ProteinData', 'cdw'))

    print("Create FBA problem")
    fba_model_file = config.get('Model', 'fba')
    biomass_reaction = config.get('Model', 'biomass_reaction')
    fba_model = create_fba_problem(fba_model_file, biomass_reaction)
    fba_model = set_bounds(fba_model, flux_data, config)
    
    print("Solve FBA problem")
    fba_model.solve()
    fluxes = {a: b for a,b in zip(fba_model.variables.get_names(),
                                  fba_model.solution.get_values())}
    biomass_reaction_flux=fluxes[biomass_reaction]

    print("Create parsimonius FBA problem")
    parsimonious_fba_model = create_parsimonious_fba_problem(sbml_file=fba_model_file, biomass_reaction=biomass_reaction,biomass_reaction_flux=biomass_reaction_flux)
    parsimonious_fba_model = set_bounds_parsimonious_fba(parsimonious_fba_model, flux_data, config)
    print("Solve parsimonious FBA problem")
    parsimonious_fba_model.solve()

    fluxes_parsimonious_fba = {a: b for a,b in zip(parsimonious_fba_model.variables.get_names(),
                                                   parsimonious_fba_model.solution.get_values())}

    net_fluxes_parsimonious_fba={}
    for rxn in fluxes_parsimonious_fba.keys():
        if rxn.endswith("__FW"):
            proto_rxn_id=rxn.split("__FW")[0]
        elif rxn.endswith("__BW"):
            proto_rxn_id=rxn.split("__BW")[0]
        else:
            proto_rxn_id=rxn

        if proto_rxn_id in fluxes_parsimonious_fba:
            net_fluxes_parsimonious_fba[proto_rxn_id]=fluxes_parsimonious_fba[proto_rxn_id]
        else:
            net_fluxes_parsimonious_fba[proto_rxn_id]=fluxes_parsimonious_fba.get("{}__FW".format(proto_rxn_id),0.0)-fluxes_parsimonious_fba.get("{}__BW".format(proto_rxn_id),0.0)

    non_zero_fluxes = dict(filter(lambda i : i[1] != 0, net_fluxes_parsimonious_fba.items()))

    print("Estimating kapp....")
    rba_model_dir = config.get('Model', 'rba')
    rba_model = rba.RbaModel.from_xml(rba_model_dir)

    gene_cnt = get_gene_cnt_per_reaction(rba_model)
    nz_gene_cnt = get_gene_cnt_per_reaction(rba_model, non_zero_fluxes)

    number_locations={p:len(protein_data.loc[p]["Cellular_Location"].split(",")) for p in protein_data.index}
    counts = {p: protein_data.loc[p][count_column]/number_locations[p] for p in protein_data.index}
    if cdw==1:
        # proteomics data are in concentration (in mmol per gcdw)
        concs = {p[0]: p[1] for p in counts.items()}
    else:
        # proteomics data are in copies per cell, and need to be converted in mmol per gcdw
        concs = {p[0]: get_mmol_gcdw(p[1], cdw) for p in counts.items()}
  

    enzymes = Enzymes(rba_model, non_zero_fluxes, counts, concs, gene_cnt, nz_gene_cnt, cdw)

    enz_ids = [enz for enz in enzymes.keys()]
    enz_kapp_df = pd.DataFrame(index=enz_ids, columns=['kapp'])
    enz_kapp_df.index.name = 'Reaction ID'

    for enz_id in enz_ids:
        enz = enzymes[enz_id]
        enz_kapp_df.loc[enz_id,'kapp']  = enz.kapp / (24*3600)

    fluxes_file = config.get('Output', 'fluxes')
    fluxes_file_type = config.get('Output', 'file_type')
    
    print("Writing fluxes....")
    flux_distribution_DF=pd.DataFrame()
    for rxn , flux in net_fluxes_parsimonious_fba.items():
    #for rxn , flux in fluxes.items():
        flux_distribution_DF.loc[rxn,"Reaction"]=rxn
        flux_distribution_DF.loc[rxn,"Flux"]=flux
    if fluxes_file_type == 'excel':
        flux_distribution_DF.to_excel(fluxes_file)
    elif fluxes_file_type == 'csv':
        delimiter = config.get('Output', 'delimiter')
        flux_distribution_DF.to_csv(fluxes_file, sep=str(delimiter))

    output_file = config.get('Output', 'kapps')
    output_file_type = config.get('Output', 'file_type')
    
    print("Writing results....")
    nz_values = enz_kapp_df
    if output_file_type == 'excel':
        nz_values.to_excel(output_file)
    elif output_file_type == 'csv':
        delimiter = config.get('Output', 'delimiter')
        nz_values.to_csv(output_file, sep=str(delimiter))
