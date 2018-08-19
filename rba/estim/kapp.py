import cobra
import ConfigParser
import pandas as pd

import rba

from utils import *


def set_bounds(fba_model, flux_data, config):
	lb_column = config.get('FluxData', 'lb_column')
	ub_column = config.get('FluxData', 'ub_column')
	for fluxes in flux_data.index:
		reactions = fluxes.split(',')
		constr = fba_model.problem.Constraint(
			sum([fba_model.reactions.get_by_id(r).flux_expression for r in reactions]),
			lb = flux_data.loc[fluxes][lb_column],
			ub = flux_data.loc[fluxes][ub_column],
		)
		fba_model.add_cons_vars(constr)


if __name__ == '__main__':
	config = ConfigParser.ConfigParser()
	config.read('kapp.cfg')

	flux_data = load_flux_data(config)

	protein_data = load_protein_data(config)
	count_column = config.get('ProteinData', 'count_column')
	cdw = float(config.get('ProteinData', 'cdw'))

	fba_model_file = config.get('Model', 'fba')
	fba_model = cobra.io.read_sbml_model(fba_model_file)

	rba_model_dir = config.get('Model', 'rba')
	rba_model = rba.RbaModel.from_xml(rba_model_dir)

	set_bounds(fba_model, flux_data, config)
	fba_sol = fba_model.optimize()
	all_fluxes = dict(fba_sol.fluxes)
	non_zero_fluxes = dict(filter(lambda i : i[1] != 0, all_fluxes.items()))

	gene_cnt = get_gene_cnt_per_reaction(rba_model)
	nz_gene_cnt = get_gene_cnt_per_reaction(rba_model, non_zero_fluxes)

	counts = {p: protein_data.loc[p][count_column] for p in protein_data.index}
	concs = {p[0]: get_mmol_gcdw(p[1], cdw) for p in counts.items()}

	enzymes = Enzymes(rba_model, non_zero_fluxes, counts, concs, gene_cnt, nz_gene_cnt, cdw)

	enz_ids = [enz for enz in enzymes.keys() if enzymes[enz].is_cytosolic]
	enz_kapp_df = pd.DataFrame(index=enz_ids, columns=['kapp'])

	for enz_id in enz_ids:
		enz = enzymes[enz_id]
		enz_kapp_df.loc[enz_id]['kapp'] = enz.kapp / 3600

	output_file = config.get('Output', 'file')
	output_file_type = config.get('Output', 'file_type')
	if output_file_type == 'excel':
		enz_kapp_df.to_excel(output_file)
	elif output_file_type == 'csv':
		delimiter = config.get('Output', 'delimiter')
		enz_kapp_df.to_csv(output_file, delimiter)
