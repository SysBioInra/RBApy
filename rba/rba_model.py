"""Module defining RbaModel class."""

# python 2/3 compatibility
from __future__ import division, print_function

# global imports
from os.path import join

# local imports
import rba.xml
from rba.core.constraint_matrix import ConstraintMatrix
from rba.core.solver import Solver
from rba.results import Results


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
        self.parameters = rba.xml.RbaParameters()
        self.proteins = rba.xml.RbaMacromolecules()
        self.enzymes = rba.xml.RbaEnzymes()
        self.rnas = rba.xml.RbaMacromolecules()
        self.dna = rba.xml.RbaMacromolecules()
        self.processes = rba.xml.RbaProcesses()
        self.medium = {}
        self.output_dir = ''

    @classmethod
    def from_xml(cls, input_dir):
        """
        Make object from xml files.

        Args:
            input_dir: path to directory containing RBA XML files.
        """
        obj = cls()
        obj.output_dir = input_dir
        obj.parameters = rba.xml.RbaParameters().from_file(
            join(input_dir, 'parameters.xml')
            )
        obj.metabolism = rba.xml.RbaMetabolism().from_file(
            join(input_dir, 'metabolism.xml')
            )
        obj.proteins = rba.xml.RbaMacromolecules().from_file(
            join(input_dir, 'proteins.xml')
            )
        obj.rnas = rba.xml.RbaMacromolecules().from_file(
            join(input_dir, 'rnas.xml')
            )
        obj.dna = rba.xml.RbaMacromolecules().from_file(
            join(input_dir, 'dna.xml')
            )
        obj.processes = rba.xml.RbaProcesses().from_file(
            join(input_dir, 'processes.xml')
            )
        obj.enzymes = rba.xml.RbaEnzymes().from_file(
            join(input_dir, 'enzymes.xml')
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

    def get_matrices(self):
        """
        Return matrix blocks used to solve the optimization problem.

        Returns:
            rba.core.RbaMatrices object containing matrix blocks.

        """
        return RbaMatrices(self)

    def solve(self, catalytic_function=None):
        """
        Solve current model.

        Args:
            catalytic_function: identifier of functions to use for catalytic
                efficiencies of enzymes ('default' if no argument given).

        Returns:
            rba.core.RbaSolver object that contains solution (if one
            was found) and matrices next to solution.

        """
        matrices = ConstraintMatrix(self)
        if catalytic_function is not None:
            matrices._blocks.set_catalytic_function(catalytic_function)
        solver = Solver(matrices)
        solver.solve()
        variables = {name: value for name, value in zip(matrices.col_names,
                                                        solver.X)}
        dual = {name: value for name, value in zip(matrices.row_names,
                                                   solver.lambda_)}
        results = Results(variables, dual, self)
        return results

    def write_files(self, output_dir=None):
        """
        Write rba files in XML format.

        Args:
            output_dir: path to directory where files should be written. If
                specified, output_dir attribute is overriden, otherwise value
                currently stored in output_dir attribute is used.
        """
        if output_dir:
            self.output_dir = output_dir
        self.metabolism.write(self._output('metabolism.xml'), 'RBAMetabolism')
        self.proteins.write(self._output('proteins.xml'), 'RBAProteins')
        self.rnas.write(self._output('rnas.xml'), 'RBARnas')
        self.dna.write(self._output('dna.xml'), 'RBADna')
        self.enzymes.write(self._output('enzymes.xml'), 'RBAEnzymes')
        self.parameters.write(self._output('parameters.xml'), 'RBAParameters')
        self.processes.write(self._output('processes.xml'), 'RBAProcesses')
        # initial conditions (medium concentrations)
        with open(self._output('medium.tsv'), 'w') as output:
            output.write('Metabolite\tConcentration\n')
            for met, conc in self.medium.items():
                output.write('{}\t{}\n'.format(met, conc))

    def _output(self, file_name):
        """Return full path to file contained in output direcotry."""
        return join(self.output_dir, file_name)
