"""Module defining RbaModel class."""

# python 2/3 compatibility
from __future__ import division, print_function, absolute_import

# global imports
from os.path import join

# local imports
import rba
from rba.utils import efficiencies


class RbaModel(object):
    """
    Class holding RBA model.

    Attributes:
        metabolism: xml structure with metabolism data.
        parameters: xml structure with parameter data.
        proteins: xml structure with protein data.
        enzymes: xml structure with enzyme data.
        rnas: xml structure with rna data.
        dna: xml structure with dna data.
        processes: xml structure with process data.
        medium: dictionary mapping metabolite prefixes with their medium
            concentration.
        output_dir: path to directory where model files should be written.

    """

    def __init__(self):
        """Build object with empty structures."""
        self.metabolism = rba.xml.RbaMetabolism()
        self.density = rba.xml.RbaDensity()
        self.parameters = rba.xml.RbaParameters()
        self.proteins = rba.xml.RbaMacromolecules()
        self.enzymes = rba.xml.RbaEnzymes()
        self.rnas = rba.xml.RbaMacromolecules()
        self.dna = rba.xml.RbaMacromolecules()
        self.processes = rba.xml.RbaProcesses()
        self.targets = rba.xml.RbaTargets()
        self.medium = {}
        self.output_dir = ''
        self._constraint_matrix = None

    @classmethod
    def from_data(cls, params_file, verbose=False):
        """
        Make object from data directory (SBML, FASTA and helper files).

        Parameters
        ----------
        params_file : str
            Path to params.in file.

        verbose : bool, optional
            Whether to display status information
        """
        builder = rba.ModelBuilder(params_file, verbose=verbose)
        builder.export_proteins('protein_summary.tsv')
        model = builder.build_model()
        return model

    @classmethod
    def from_xml(cls, input_dir):
        """
        Make object from xml files.

        Parameters
        ----------
        input_dir : str
            Path to directory containing RBA XML files.

        """
        obj = cls()
        obj.output_dir = input_dir
        obj.parameters = rba.xml.RbaParameters().from_file(
            open(join(input_dir, 'parameters.xml'))
        )
        obj.density = rba.xml.RbaDensity().from_file(
            open(join(input_dir, 'density.xml'))
        )
        obj.metabolism = rba.xml.RbaMetabolism().from_file(
            open(join(input_dir, 'metabolism.xml'))
        )
        obj.proteins = rba.xml.RbaMacromolecules().from_file(
            open(join(input_dir, 'proteins.xml'))
        )
        obj.rnas = rba.xml.RbaMacromolecules().from_file(
            open(join(input_dir, 'rnas.xml'))
        )
        obj.dna = rba.xml.RbaMacromolecules().from_file(
            open(join(input_dir, 'dna.xml'))
        )
        obj.processes = rba.xml.RbaProcesses().from_file(
            open(join(input_dir, 'processes.xml'))
        )
        obj.targets = rba.xml.RbaTargets().from_file(
            open(join(input_dir, 'targets.xml'))
        )
        obj.enzymes = rba.xml.RbaEnzymes().from_file(
            open(join(input_dir, 'enzymes.xml'))
        )
        obj.set_medium(join(input_dir, 'medium.tsv'))
        return obj

    def set_medium(self, file_name):
        """
        Set medium concentrations according to file.

        Args:
            file_name: path to file containing medium concentrations in a tab
                separated format. File is supposed to contain a header, one
                column with metabolite prefixes and one column with
                concentrations values.
        """
        concentrations = {}
        with open(file_name, 'rU') as input_stream:
            # skip header
            next(input_stream)
            for line in input_stream:
                [met, conc] = line.rstrip().split('\t')
                concentrations[met] = float(conc)
        self.medium = concentrations
        if self._constraint_matrix:
            self._constraint_matrix.set_medium(concentrations)

    def write(self, output_dir=None):
        """
        Write rba files in XML format.

        Parameters
        ----------
        output_dir : str
            Path to directory where files should be written. If
            specified, output_dir attribute is overriden, otherwise value
            currently stored in output_dir attribute is used.

        """
        if output_dir:
            self.output_dir = output_dir
        self.metabolism.write(self._output('metabolism.xml'), 'RBAMetabolism')
        self.density.write(self._output('density.xml'), 'RBADensity')
        self.proteins.write(self._output('proteins.xml'), 'RBAProteins')
        self.rnas.write(self._output('rnas.xml'), 'RBARnas')
        self.dna.write(self._output('dna.xml'), 'RBADna')
        self.enzymes.write(self._output('enzymes.xml'), 'RBAEnzymes')
        self.parameters.write(self._output('parameters.xml'), 'RBAParameters')
        self.processes.write(self._output('processes.xml'), 'RBAProcesses')
        self.targets.write(self._output('targets.xml'), 'RBATargets')
        # initial conditions (medium concentrations)
        with open(self._output('medium.tsv'), 'w') as output:
            output.write('Metabolite\tConcentration\n')
            for met, conc in self.medium.items():
                output.write('{}\t{}\n'.format(met, conc))

    def _output(self, file_name):
        """Return full path to file contained in output directory."""
        return join(self.output_dir, file_name)

    def solve(self, recompute_matrices=True, lp_solver=None,
              mu_min=0., mu_max=2.5,
              bissection_tol=1e-6, max_bissection_iters=None,
              verbose=False):
        """
        Solve RBA model.

        Parameters
        ----------
        recompute_matrices : bool, optional
            If the model is solved several time, recompute matrices defining
            the optimality problem (True by default). This parameter should be
            set to False when only the medium composition changes (medium
            concentrations do not appear in the matrices).

        lp_solver : str, optional
            LP solver (``cplex``, ``glpk``, ``gurobi``, ``scipy`` ``swiglpk``)

        mu_min : float, optional:
            Minimum μ to check

        mu_max : float, optional:
            Maximum μ to check

        bissection_tol : float, optional
            Tolerance for bissection

        max_bissection_iters : int, optional
            Maximum number of iterations for bissection

        verbose : bool, optional
            Whether to display status information

        Returns
        -------
        rba.utils.Results object containing optimal growth rate and fluxes.

        """
        if recompute_matrices or self._constraint_matrix is None:
            self._constraint_matrix = rba.ConstraintMatrix(self)
        solver = rba.Solver(self._constraint_matrix,
                            lp_solver=lp_solver,
                            mu_min=mu_min,
                            mu_max=mu_max,
                            bissection_tol=bissection_tol,
                            max_bissection_iters=max_bissection_iters,
                            verbose=verbose)
        solver.solve()
        # solver.solve_grid()
        return rba.Results(self, self._constraint_matrix, solver)

    def set_enzyme_efficiencies(self, file_name):
        efficiencies.set_efficiencies(self, file_name)
