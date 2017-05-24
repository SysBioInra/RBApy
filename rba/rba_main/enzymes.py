
import numpy

from .functions import *

class Enzymes(object):
    def __init__(self, enzymes, species, reactions):
        # check that all reactions are found and keep only enzymes
        # that have a machinery
        reactions_left = reactions[:]
        nonzero_enzymes = []
        self.reaction_catalyzed = []
        for e in enzymes.enzymes:
            reactions_left.remove(e.enzymatic_activity.reaction)
            if not(e.zero_cost or e.machinery_composition.is_empty()):
                nonzero_enzymes.append(e)
                self.reaction_catalyzed.append(e.enzymatic_activity.reaction)
        self.ids = [e.id for e in nonzero_enzymes]        
        if len(reactions_left) > 0:
            print('Warning: did not find enzymes for following reactions: '
                  + ', '.join(reactions_left))
            
        # extract machinery information
        machinery = [e.machinery_composition for e in nonzero_enzymes]
        self.machinery = species.create_machinery(machinery)
        
        # extract efficiency information
        self.efficiency = EnzymeEfficiency(enzymes.efficiency_functions)
        for enzyme in nonzero_enzymes:
            self.efficiency.add_activity(enzyme.enzymatic_activity)

class EnzymeEfficiency(object):
    def __init__(self, efficiency_functions):
        # read list of efficiency functions
        self._fn_types = {}
        for fn in efficiency_functions: self._fn_types[fn.id] = fn.type
        self._efficiency = {fn_id: [] for fn_id in self._fn_types}
        self._import = []
        # use default values for import and efficiency function
        self._eff_fn = self._efficiency[efficiency_functions[0].id]
        self._import_values = 1

    def add_activity(self, enzymatic_activity):
        # base efficiency
        params = {}
        for fn in enzymatic_activity.enzyme_efficiencies:
            params[fn.function] = {p.id: p.value for p in fn.parameters}
        for fn_id in self._efficiency:
            try:
                self._efficiency[fn_id].append \
                    (build_function(self._fn_types[fn_id], params[fn_id]))
            except KeyError:
                print('Missing parameters for enzymatic function ' + fn_id)
                raise UserWarning('Invalid enzyme file.')
        # import efficiency
        t_eff = enzymatic_activity.transporter_efficiency
        if t_eff:
            fns = []
            for fn in t_eff:
                params = {p.id: p.value for p in fn.parameters}
                fns.append(build_function(fn.type, params, fn.variable))
            self._import.append(fns)
        else:
            self._import.append([])

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
