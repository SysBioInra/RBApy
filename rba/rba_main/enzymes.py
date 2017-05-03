
import numpy

from functions import *

class Enzymes(object):
    def __init__(self, data, species):
        # extract machinery information
        machinery = [e.machinery_composition for e in data.enzymes.enzymes]
        self.machinery = species.create_machinery(machinery)
        # extract efficiency information
        self.efficiency = EnzymeEfficiency(data)

    def add_null_enzymes(self, number):
        self.efficiency.add_null_enzymes(number)
        self.machinery.add_null_machineries(number)

class EnzymeEfficiency(object):
    def __init__(self, data):
        # read list of efficiency functions
        fn_types = {}
        for fn in data.enzymes.efficiency_functions: fn_types[fn.id] = fn.type
        # read efficiencies
        self._efficiency = {fn_id: [] for fn_id in fn_types}
        self._import = []
        for enzyme in data.enzymes.enzymes:
            # base efficiency
            params = {}
            for fn in enzyme.enzymatic_activity.enzyme_efficiencies:
                params[fn.function] = {p.id: p.value for p in fn.parameters}
            for fn_id in self._efficiency:
                try:
                    self._efficiency[fn_id].append \
                        (build_function(fn_types[fn_id], params[fn_id]))
                except KeyError:
                    print('Missing parameters for enzymatic function ' + fn_id)
                    raise UserWarning('Invalid enzyme file.')
            # import efficiency
            t_eff = enzyme.enzymatic_activity.transporter_efficiency
            if t_eff:
                fns = []
                for fn in t_eff:
                    params = {p.id: p.value for p in fn.parameters}
                    fns.append(build_function(fn.type, params, fn.variable))
                self._import.append(fns)
            else:
                self._import.append([])
        # use default values for import and efficiency function
        self._eff_fn = self._efficiency[data.enzymes.efficiency_functions[0].id]
        self._import_values = numpy.ones(len(self._import))

    def add_null_enzymes(self, number):
        for eff_fn in self._efficiency.itervalues(): eff_fn += [None]*number
        self._import += [[]]*number
        self._import_values = numpy.append(self._import_values,
                                           numpy.ones(number))

    def set_function(self, fn_id):
        self._eff_fn = self._efficiency[fn_id]
        
    def update_import(self, concentration):
        self._import_values = numpy.ones(len(self._import))
        for i, import_fn in enumerate(self._import):
            for i_fn in import_fn:
                # /!\ we identify metabolites by their prefix !!!
                key = i_fn.variable.rsplit('_',1)[0]
                self._import_values[i] *= i_fn.evaluate(concentration[key])

    def compute(self, mu):
        efficiency = numpy.ones(len(self._eff_fn))
        for i, fn in enumerate(self._eff_fn):
            if fn: efficiency[i] = fn.evaluate(mu)
        return (efficiency*self._import_values, efficiency)
