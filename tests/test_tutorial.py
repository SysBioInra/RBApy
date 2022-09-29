""" Tutorial (``docs/xml_tutorial.ipynb``) encoded into unit tests
"""

from rba.core import solver
import numpy.testing
import os
import rba
import shutil
import tempfile
import unittest


class TutorialTestCase(unittest.TestCase):
    def setUp(self):
        self.tmp_dirname = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmp_dirname)

    @unittest.skipIf(not solver.is_cplex_available(), reason='CPLEX is not available')
    def test_cplex(self):
        self._test('cplex', rtol=1e-7)

    # @unittest.skipIf(
    #    not (solver.is_conv_opt_available() and solver.is_cplex_available()),
    #    reason='ConvOpt and/or CPLEX are not available')
    # def test_cplex_conv_opt(self):
    #    self._test('cplex_conv_opt', rtol=1e-7)

    # @unittest.skipIf(
    #    not (solver.is_optlang_available() and solver.is_cplex_available()),
    #    reason='OptLang and/or CPLEX are not available')
    # def test_cplex_optlang(self):
    #    self._test('cplex_optlang', rtol=1e-7)

    @unittest.skipIf(
        not (solver.is_glpk_available()),
        reason='GLPK is not available')
    def test_glpk(self):
        self._test('glpk')

    @unittest.skipIf(
        not (solver.is_swiglpk_available()),
        reason='SWIGLPK is not available')
    def test_swiglpk(self):
        self._test('swiglpk')

    @unittest.skipIf(
        not (solver.is_optlang_available() and solver.is_swiglpk_available()),
        reason='OptLang and/or GLPK are not available')
    def test_glpk_optlang(self):
        self._test('glpk_optlang')

    @unittest.skipIf(
        not (solver.is_gurobi_available()),
        reason='Gurobi is not available')
    def test_gurobi(self):
        self._test('gurobi')

    @unittest.skipIf(
        not (solver.is_conv_opt_available() and solver.is_gurobi_available()),
        reason='ConvOpt and/or Gurobi are not available')
    def test_gurobi_conv_opt(self):
        self._test('gurobi_conv_opt')

    # @unittest.skipIf(
    #     not (solver.is_optlang_available() and solver.is_gurobi_available()),
    #     reason='OptLang and/or Gurobi are not available')
    # def test_gurobi_optlang(self):
    #    self._test('gurobi_optlang')

    def _test(self, lp_solver, rtol=5e-3):
        ######################################################
        # Creation of an empty model with RBApy
        my_model = rba.RbaModel()
        self.assertEqual(len(my_model.metabolism.compartments), 0)

        my_model.write(self.tmp_dirname)

        my_model = rba.RbaModel.from_xml(self.tmp_dirname)
        self.assertEqual(len(my_model.metabolism.compartments), 0)

        ######################################################
        # metabolism.xml: compartments, metabolites and reactions
        my_model.metabolism.species.append(rba.xml.Species('M_carbon_source_e', boundary_condition=True))
        my_model.metabolism.species.append(rba.xml.Species('M_carbon_source_c', boundary_condition=False))
        my_model.metabolism.species.append(rba.xml.Species('M_protein_component_precursor_c', boundary_condition=False))
        my_model.metabolism.species.append(rba.xml.Species('M_biomass_c', boundary_condition=False))

        self.assertEqual(len(my_model.metabolism.species), 4)
        self.assertEqual(my_model.metabolism.species[0].id, 'M_carbon_source_e')
        # my_model.write()

        my_model.metabolism.compartments.append(rba.xml.Compartment('extracellular'))
        my_model.metabolism.compartments.append(rba.xml.Compartment('cytosol'))
        # my_model.write()

        reaction_1 = rba.xml.Reaction('R_transport', reversible=False)
        reaction_1.reactants.append(rba.xml.SpeciesReference('M_carbon_source_e', 1))
        reaction_1.products.append(rba.xml.SpeciesReference('M_carbon_source_c', 1))
        reaction_2 = rba.xml.Reaction('R_protein_component_precursor', reversible=False)
        reaction_2.reactants.append(rba.xml.SpeciesReference('M_carbon_source_c', 1))
        reaction_2.products.append(rba.xml.SpeciesReference('M_protein_component_precursor_c', 1))
        reaction_3 = rba.xml.Reaction('R_biomass', reversible=False)
        reaction_3.reactants.append(rba.xml.SpeciesReference('M_carbon_source_c', 1))
        reaction_3.products.append(rba.xml.SpeciesReference('M_biomass_c', 1))
        my_model.metabolism.reactions.append(reaction_1)
        my_model.metabolism.reactions.append(reaction_2)
        my_model.metabolism.reactions.append(reaction_3)
        # my_model.write()

        my_model.medium['M_carbon_source'] = 10
        # my_model.write()

        results = my_model.solve(lp_solver=lp_solver)
        numpy.testing.assert_allclose(results.mu_opt, 2.5, rtol=rtol)
        self.assertEqual(results.variables, {
            'R_transport': -0.0,
            'R_protein_component_precursor': 0.0,
            'R_biomass': 0.0,
        })
        self.assertEqual(results.dual_values, {
            'M_carbon_source_c': 0.0,
            'M_protein_component_precursor_c': 0.0,
            'M_biomass_c': 0.0,
        })

        ######################################################
        # proteins.xml, rnas.xml, dna.xml: macromolecules
        my_model.proteins.components.append(rba.xml.Component(id_='protein_component_residue',
                                                              name='Protein residue', type_='residue', weight=1))

        protein_1 = rba.xml.Macromolecule('small_protein', 'cytosol')
        protein_1.composition.append(rba.xml.ComponentReference('protein_component_residue', 10))
        my_model.proteins.macromolecules.append(protein_1)
        protein_2 = rba.xml.Macromolecule('large_protein', 'cytosol')
        protein_2.composition.append(rba.xml.ComponentReference('protein_component_residue', 20))
        my_model.proteins.macromolecules.append(protein_2)

        # my_model.write()

        ######################################################
        # enzymes.xml: catalytic constraints
        enzyme_1 = rba.xml.Enzyme(id_='R_transport_enzyme', reaction='R_transport',
                                  forward_efficiency='kcat_transport', backward_efficiency='zero')
        enzyme_1.machinery_composition.reactants.append(rba.xml.SpeciesReference('large_protein', 2))
        my_model.enzymes.enzymes.append(enzyme_1)

        enzyme_2 = rba.xml.Enzyme(id_='R_protein_precursor_enzyme', reaction='R_protein_component_precursor',
                                  forward_efficiency='kcat_precursor', backward_efficiency='zero')
        enzyme_2.machinery_composition.reactants.append(rba.xml.SpeciesReference('small_protein', 2))
        my_model.enzymes.enzymes.append(enzyme_2)
        enzyme_3 = rba.xml.Enzyme(id_='R_biomass_enzyme', reaction='R_biomass',
                                  forward_efficiency='kcat_biomass', backward_efficiency='zero')
        enzyme_3.machinery_composition.reactants.append(rba.xml.SpeciesReference('small_protein', 1))
        enzyme_3.machinery_composition.reactants.append(rba.xml.SpeciesReference('large_protein', 1))
        my_model.enzymes.enzymes.append(enzyme_3)

        # my_model.write()

        ######################################################
        # parameters.xml: centralizing all model parameters
        parameter_1 = rba.xml.Function('kcat_transport_base', 'constant', {'CONSTANT': 10})
        my_model.parameters.functions.append(parameter_1)
        parameter_2 = rba.xml.Function('kcat_precursor', 'constant', {'CONSTANT': 10})
        my_model.parameters.functions.append(parameter_2)
        parameter_3 = rba.xml.Function('kcat_biomass', 'constant', {'CONSTANT': 10})
        my_model.parameters.functions.append(parameter_3)
        parameter_4 = rba.xml.Function('zero', 'constant', {'CONSTANT': 0})
        my_model.parameters.functions.append(parameter_4)

        parameter_4 = rba.xml.Function('transport_factor', 'michaelisMenten', {'kmax': 1, 'Km': 0.5}, variable='M_carbon_source_e')
        my_model.parameters.functions.append(parameter_4)

        aggregate_1 = rba.xml.Aggregate('kcat_transport', type_='multiplication')
        aggregate_1.function_references.append(rba.xml.FunctionReference('kcat_transport_base'))
        aggregate_1.function_references.append(rba.xml.FunctionReference('transport_factor'))
        my_model.parameters.aggregates.append(aggregate_1)

        # my_model.write()

        results = my_model.solve(lp_solver=lp_solver)
        numpy.testing.assert_allclose(results.mu_opt, 2.5, rtol=rtol)
        """
        self.assertEqual(results.variables, {
            'R_transport': 0.0,
            'R_protein_component_precursor': 0.0,
            'R_biomass': 0.0,
            'R_transport_enzyme': -0.0,
            'R_protein_precursor_enzyme': 0.0,
            'R_biomass_enzyme': 0.0
        })
        self.assertEqual(results.dual_values, {
            'M_carbon_source_c': 0.0,
            'M_protein_component_precursor_c': 0.0,
            'M_biomass_c': 0.0,
            'R_transport_enzyme_forward_capacity': -0.105,
            'R_protein_precursor_enzyme_forward_capacity': 0.0,
            'R_biomass_enzyme_forward_capacity': 0.0,
            'R_transport_enzyme_backward_capacity': 0.0,
            'R_protein_precursor_enzyme_backward_capacity': 0.0,
            'R_biomass_enzyme_backward_capacity': 0.0,
        })
        """

        ######################################################
        # density.xml: density constraints for compartments

        density_constraint = rba.xml.TargetDensity('cytosol')
        density_constraint.upper_bound = 'maximal_cytosol_density'
        my_model.density.target_densities.append(density_constraint)

        my_model.parameters.functions.append(rba.xml.Function('maximal_cytosol_density', 'constant', {'CONSTANT': 10}))

        # my_model.write()

        results = my_model.solve(lp_solver=lp_solver)
        numpy.testing.assert_allclose(results.mu_opt, 2.5, rtol=rtol)
        """
        self.assertEqual(results.variables, {
            'R_transport': 0.0,
            'R_protein_component_precursor': 0.0,
            'R_biomass': 0.0,
            'R_transport_enzyme': -0.0,
            'R_protein_precursor_enzyme': 0.0,
            'R_biomass_enzyme': 0.0
        })
        self.assertEqual(results.dual_values, {
            'M_carbon_source_c': 0.0,
            'M_protein_component_precursor_c': 0.0,
            'M_biomass_c': 0.0,
            'R_transport_enzyme_forward_capacity': -0.105,
            'R_protein_precursor_enzyme_forward_capacity': 0.0,
            'R_biomass_enzyme_forward_capacity': 0.0,
            'R_transport_enzyme_backward_capacity': 0.0,
            'R_protein_precursor_enzyme_backward_capacity': 0.0,
            'R_biomass_enzyme_backward_capacity': 0.0,
            'cytosol_density': 0.0
        })
        """

        ######################################################
        # targets.xml: production requirements
        new_target = rba.xml.TargetSpecies('M_biomass_c')
        new_target.value = 'target_biomass_production'
        new_target_group = rba.xml.TargetGroup('biomass_production')
        new_target_group.concentrations.append(new_target)
        my_model.targets.target_groups.append(new_target_group)

        my_model.parameters.functions.append(rba.xml.Function('target_biomass_production', 'constant', {'CONSTANT': 1}))

        # my_model.write()

        results = my_model.solve(lp_solver=lp_solver)
        numpy.testing.assert_allclose(results.mu_opt, 1.3888883590698242, rtol=rtol)
        """
        self.assertEqual(results.variables, {
            'R_transport': 1.3888883590698242,
            'R_protein_component_precursor': 0.0,
            'R_biomass': 1.3888883590698242,
            'R_transport_enzyme': 0.14583327770233154,
            'R_protein_precursor_enzyme': 0.0,
            'R_biomass_enzyme': 0.13888883590698242
        })
        self.assertEqual(results.dual_values, {
            'M_carbon_source_c': 0.105,
            'M_protein_component_precursor_c': 0.0,
            'M_biomass_c': 0.20500000000000002,
            'R_transport_enzyme_forward_capacity': -0.105,
            'R_protein_precursor_enzyme_forward_capacity': 0.0,
            'R_biomass_enzyme_forward_capacity': -0.1,
            'R_transport_enzyme_backward_capacity': 0.0,
            'R_protein_precursor_enzyme_backward_capacity': 0.0,
            'R_biomass_enzyme_backward_capacity': 0.0,
            'cytosol_density': 0.0
        })
        """

        maximal_density_fn = my_model.parameters.functions.get_by_id('maximal_cytosol_density')
        maximal_density = maximal_density_fn.parameters.get_by_id('CONSTANT')
        maximal_density.value = 10
        results = my_model.solve(lp_solver=lp_solver)
        numpy.testing.assert_allclose(results.mu_opt, 1.3888883590698242, rtol=rtol)
        maximal_density.value = 1
        results = my_model.solve(lp_solver=lp_solver)
        numpy.testing.assert_allclose(results.mu_opt, 0.13888835906982422, rtol=rtol)

        ######################################################
        # processes.xml: macromolecule production
        new_map = rba.xml.ProcessingMap(id_='translation_map')
        new_map.constant_processing.reactants.append(rba.xml.SpeciesReference('M_carbon_source_c', 1))
        residue_processing = rba.xml.ComponentProcessing(component='protein_component_residue')
        residue_processing.reactants.append(rba.xml.SpeciesReference('M_protein_component_precursor_c', 1))
        residue_processing.reactants.append(rba.xml.SpeciesReference('M_carbon_source_c', 1))
        new_map.component_processings.append(residue_processing)
        my_model.processes.processing_maps.append(new_map)

        translation = rba.xml.Process(id_='translation', name='Translation process')
        processing = rba.xml.Processing(map_='translation_map', set_='protein')
        for protein in my_model.proteins.macromolecules:
            processing.inputs.append(rba.xml.SpeciesReference(protein.id, 1))
        translation.processings.productions.append(processing)
        my_model.processes.processes.append(translation)

        # my_model.write()

        results = my_model.solve(lp_solver=lp_solver)
        numpy.testing.assert_allclose(results.mu_opt, 0.05603671073913574, rtol=rtol)

        ribosome = rba.xml.Machinery()
        ribosome.machinery_composition.reactants.append(rba.xml.SpeciesReference('small_protein', 1))
        ribosome.machinery_composition.reactants.append(rba.xml.SpeciesReference('large_protein', 1))
        my_model.processes.processes.get_by_id('translation').machinery = ribosome

        translation_map = my_model.processes.processing_maps.get_by_id('translation_map')
        translation_map.component_processings[0].machinery_cost = 0
        ribosome.capacity.value = 'ribosome_capacity'
        my_model.parameters.functions.append(rba.xml.Function('ribosome_capacity', 'constant', {'CONSTANT': 100}))

        my_model.write()

        results = my_model.solve(lp_solver=lp_solver)
        numpy.testing.assert_allclose(results.mu_opt, 0.05603671073913574, rtol=rtol)

        ######################################################
        # Parameters.xml: overview of the impact of parameters
        maximal_density = my_model.parameters.functions.get_by_id('maximal_cytosol_density').parameters.get_by_id('CONSTANT')
        transport_kcat = my_model.parameters.functions.get_by_id('kcat_transport_base').parameters.get_by_id('CONSTANT')
        precursor_kcat = my_model.parameters.functions.get_by_id('kcat_precursor').parameters.get_by_id('CONSTANT')
        biomass_kcat = my_model.parameters.functions.get_by_id('kcat_biomass').parameters.get_by_id('CONSTANT')
        ribosome_capacity = my_model.parameters.functions.get_by_id('ribosome_capacity').parameters.get_by_id('CONSTANT')

        maximal_density.value = 10
        results = my_model.solve(lp_solver=lp_solver)
        numpy.testing.assert_allclose(results.mu_opt, 0.08795976638793945, rtol=rtol)
        maximal_density.value = 100
        results = my_model.solve(lp_solver=lp_solver)
        numpy.testing.assert_allclose(results.mu_opt, 0.09327113628387451, rtol=rtol)
        maximal_density.value = 1000
        results = my_model.solve(lp_solver=lp_solver)
        numpy.testing.assert_allclose(results.mu_opt, 0.09383797645568848, rtol=rtol)
        maximal_density.value = 10000
        results = my_model.solve(lp_solver=lp_solver)
        numpy.testing.assert_allclose(results.mu_opt, 0.09389519691467285, rtol=rtol)

        transport_kcat.value = 100
        precursor_kcat.value = 100
        biomass_kcat.value = 100
        results = my_model.solve(lp_solver=lp_solver)
        numpy.testing.assert_allclose(results.mu_opt, 0.938953161239624, rtol=rtol)

        precursor_kcat.value = 100000
        biomass_kcat.value = 100000
        results = my_model.solve(lp_solver=lp_solver)
        numpy.testing.assert_allclose(results.mu_opt, 1.161106824874878, rtol=rtol)
        precursor_kcat.value = 1000000
        biomass_kcat.value = 1000000
        results = my_model.solve(lp_solver=lp_solver)
        numpy.testing.assert_allclose(results.mu_opt, 1.1613553762435913, rtol=rtol)

        transport_kcat.value = 100
        results = my_model.solve(lp_solver=lp_solver)
        numpy.testing.assert_allclose(results.mu_opt, 1.1613553762435913, rtol=rtol)
        transport_kcat.value = 200
        results = my_model.solve(lp_solver=lp_solver)
        numpy.testing.assert_allclose(results.mu_opt, 2.322656512260437, rtol=rtol)
        transport_kcat.value = 50
        results = my_model.solve(lp_solver=lp_solver)
        numpy.testing.assert_allclose(results.mu_opt, 0.5806845426559448, rtol=rtol)

        maximal_density.value = 10000
        results = my_model.solve(lp_solver=lp_solver)
        numpy.testing.assert_allclose(results.mu_opt, 0.5806845426559448, rtol=rtol)
        maximal_density.value = 100
        results = my_model.solve(lp_solver=lp_solver)
        numpy.testing.assert_allclose(results.mu_opt, 0.577893853187561, rtol=rtol)
        maximal_density.value = 10
        results = my_model.solve(lp_solver=lp_solver)
        numpy.testing.assert_allclose(results.mu_opt, 0.5537021160125732, rtol=rtol)
        maximal_density.value = 1
        results = my_model.solve(lp_solver=lp_solver)
        numpy.testing.assert_allclose(results.mu_opt, 0.3903120756149292, rtol=rtol)
