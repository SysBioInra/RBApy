from rba.core import solver
import rba.cli.generate_rba_model
import rba.cli.solve_rba_model
import os
import rba
import shutil
import tempfile
import unittest
import unittest.mock


class CliTestCase(unittest.TestCase):
    def setUp(self):
        self.tmp_dirname = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmp_dirname)

    def test(self):
        # help
        with unittest.mock.patch('sys.argv', ['', '--help']):
            with self.assertRaises(SystemExit):
                rba.cli.generate_rba_model.main()

        with unittest.mock.patch('sys.argv', ['', '--help']):
            with self.assertRaises(SystemExit):
                rba.cli.solve_rba_model.main()

        # model generation
        parameters_file = os.path.join(os.path.dirname(__file__), '..', 'sample_input_small', 'params.in')
        with unittest.mock.patch('sys.argv', ['', parameters_file, '--model-dir', self.tmp_dirname]):
            rba.cli.generate_rba_model.main()

        self.assertTrue(os.path.join(self.tmp_dirname, 'metabolism.xml'))

        parameters_file = os.path.join(os.path.dirname(__file__), '..', 'sample_input_small', 'params.in')
        with unittest.mock.patch('sys.argv', ['', parameters_file, '--model-dir', self.tmp_dirname, '--verbose']):
            rba.cli.generate_rba_model.main()

        # simulation
        if (
            solver.is_cplex_available()
            or solver.is_glpk_available()
            or solver.is_swiglpk_available()
            or solver.is_gurobi_available()
        ):
            with unittest.mock.patch('sys.argv', ['', self.tmp_dirname, '--bissection-tol', '1e-2', '--verbose']):
                results = rba.cli.solve_rba_model.main()

            with unittest.mock.patch('sys.argv', ['', self.tmp_dirname, '--bissection-tol', '1e-2']):
                results = rba.cli.solve_rba_model.main()

        """
        with unittest.mock.patch('sys.argv', ['',
                                              self.tmp_dirname,
                                              '--lp-solver', 'glpk',
                                              '--bissection-tol', '1e-2',
                                              '--output-dir', os.path.join(self.tmp_dirname, 'results'),
                                              '--verbose',
                                              ]):
            results = rba.cli.solve_rba_model.main()
        """
