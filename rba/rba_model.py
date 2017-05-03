
import os.path

from rba_xml import *
from rba_main.rba_matrices import RbaMatrices
from rba_main.rba_solver import RbaSolver

class RbaModel(object):
    def __init__(self):
        self.metabolism = RbaMetabolism()
        self.parameters = RbaParameters()
        self.proteins = RbaMacromolecules()
        self.enzymes = RbaEnzymes()
        self.rnas = RbaMacromolecules()
        self.dna = RbaMacromolecules()
        self.processes = RbaProcesses()
        self.medium = {}
        self.output_dir = ''

    @classmethod
    def from_xml(cls, input_dir):
        obj = cls()
        obj._input_dir = input_dir
        obj.output_dir = input_dir
        obj.parameters = RbaParameters().from_file(obj._input('parameters.xml'))
        obj.metabolism = RbaMetabolism().from_file(obj._input('metabolism.xml'))
        obj.proteins = RbaMacromolecules().from_file(obj._input('proteins.xml'))
        obj.rnas = RbaMacromolecules().from_file(obj._input('rnas.xml'))
        obj.dna = RbaMacromolecules().from_file(obj._input('dna.xml'))
        obj.processes = RbaProcesses().from_file(obj._input('processes.xml'))
        obj.enzymes = RbaEnzymes().from_file(obj._input('enzymes.xml'))
        obj.medium = obj.read_medium(obj._input('medium.tsv'))
        return obj
                                       
    def read_medium(self, file_name):
        concentrations = {}
        with open(file_name, 'r') as f:
            # skip header
            next(f)
            for line in f:
                [met, conc] = line.rstrip().split('\t')
                concentrations[met] = float(conc)
        return concentrations

    def solve(self, catalytic_function = None, scaling_factor = 1000):
        matrices = RbaMatrices(self)
        if catalytic_function is not None:
            matrices.enzymes.efficiency.set_function(catalytic_function)
        return RbaSolver(matrices).solve(scaling_factor)

    def write_files(self, output_dir = None):
        """
        Write rba files.
        """
        if output_dir: self.output_dir = output_dir
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
            for m, conc in self.medium.iteritems():
                output.write(m + '\t' + str(conc) + '\n')

    def _input(self, file_name):
        return os.path.join(self._input_dir, file_name)
        
    def _output(self, file_name):
        return os.path.join(self.output_dir, file_name)
