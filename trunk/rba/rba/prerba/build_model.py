"""Module with main function importing data and building model."""

# python 2/3 compatibility
from __future__ import division, print_function, absolute_import

# global imports
import os.path
import pandas

# local import
from rba.prerba.pipeline_parameters import PipelineParameters
from rba.rba_model import RbaModel
from rba.prerba.model_builder import ModelBuilder
from rba.prerba.default_data import DefaultData
from rba.prerba.user_data import UserData
from rba.prerba import protein_export


def build_model(parameter_file):
    """
    Generate RBA model from SBML file and helper files.

    Parameters
    ----------
    parameter_file: str
        standard pipeline input file (containing e.g.
        sbml location, organism id).

    """
    parameters = PipelineParameters(parameter_file).parameters
    default_data = DefaultData()
    user_data = UserData(parameters, default_data)
    protein_export.export_proteins(
        os.path.join(parameters['INPUT_DIR'], 'protein_summary.tsv'),
        user_data.enzymatic_proteins
        )
    print('Building model...')
    builder = ModelBuilder(default_data, user_data)
    model = RbaModel()
    model.metabolism = builder.build_metabolism()
    model.density = builder.build_density()
    model.parameters = builder.build_parameters()
    model.proteins = builder.build_proteins()
    model.rnas = builder.build_rnas()
    model.dna = builder.build_dna()
    model.processes = builder.build_processes()
    model.enzymes = builder.build_enzymes()
    concentration = default_data.activity.medium_concentration
    for metab in user_data.external_metabolites():
        # !!! we identify metabolites by their prefix !!!
        # the idea is that M_glc_p and M_glc_e will be seen
        # as the same metabolite M_glc
        model.medium[metab.rsplit('_', 1)[0]] = concentration
    model.output_dir = parameters['OUTPUT_DIR']
    print('Done.')
    return model
