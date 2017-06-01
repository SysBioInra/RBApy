
from scipy.sparse import hstack, csr_matrix, lil_matrix
from itertools import chain

from .target_vector import TargetVector

class Processes(object):
    def __init__(self, processes, species, functions):
        ## extract ids
        self.ids = [p.id for p in processes]

        ## extract machinery related information
        # extract machinery
        machinery = [p.machinery.machinery_composition for p in processes]
        self.machinery = species.create_machinery(machinery)
        # extract capacities
        # we use CPLEX conventions for signs E=equality, L=lower than.
        values = []
        self.capacity_signs = []
        for p in processes:
            if p.machinery.capacity.value is not None:
                values.append(p.machinery.capacity.value)
                self.capacity_signs.append('E')
            elif p.machinery.capacity.upper_bound is not None:
                values.append(p.machinery.capacity.value)
                self.capacity_signs.append('L')
            else:
                values.append(0)
                self.capacity_signs.append('E')
        self.capacity = TargetVector(values, functions)
        
        ## extract targets
        # extract target values (absolute + concentration related fluxes)
        self.target_values = TargetValues(processes, species, functions)
        self.undetermined_values = UndeterminedValues(processes, species,
                                                      functions)
        # extract target reactions
        self.target_reaction = []
        self.lb_reaction = []
        self.ub_reaction = []
        val = []
        lb = []
        ub = []
        for p in processes:
            for r in p.targets.reaction_fluxes:
                if r.value is not None:
                    val.append(r.value)
                    self.target_reaction.append(r.reaction)
                else:
                    if r.lower_bound is not None:
                        lb.append(r.lower_bound)
                        self.lb_reaction.append(r.reaction)
                    if r.upper_bound is not None:
                        ub.append(r.upper_bound)
                        self.ub_reaction.append(r.reaction)
        self.reaction_value = TargetVector(val, functions)
        self.lb = TargetVector(lb, functions)
        self.ub = TargetVector(ub, functions)

class TargetValues(object):
    def __init__(self, processes, species, functions):
        targets = [p.targets for p in processes]
        # read fluxes used to maintain concentration
        # (only if a value is specified!)
        target_species = [s for t in targets for s in t.concentrations
                          if s.value is not None]
        indices = [species.ids.index(t.species) for t in target_species];
        self._conc_values = TargetVector([s.value for s in target_species],
                                         functions)
        self._conc_composition = species.production[:, indices].tocsr()
        self._conc_proc_cost = species.prod_proc_cost[:, indices].tocsr()
        self._weight = species.weight[:, indices].tocsr()
        # read production and degradation fluxes
        # (only if a value is specified!)
        prod_targets = [s for t in targets for s in t.production_fluxes
                        if s.value is not None]
        prod_indices = [species.ids.index(t.species) for t in prod_targets];
        deg_targets = [s for t in targets for s in t.degradation_fluxes
                       if s.value is not None]
        deg_indices = [species.ids.index(t.species) for t in deg_targets];
        self._values = TargetVector([s.value
                                     for s in chain(prod_targets,deg_targets)],
                                    functions)
        self._composition = hstack([species.production[:,prod_indices],
                                    species.degradation[:,deg_indices]]).tocsr()
        self._proc_cost = hstack([species.prod_proc_cost[:,prod_indices],
                                  species.deg_proc_cost[:,deg_indices]]).tocsr()

    def compute(self, mu):
        conc_values = self._conc_values.compute(mu)
        abs_values = self._values.compute(mu)
        return (self._conc_composition * (mu * conc_values)
                + self._composition * abs_values,
                self._conc_proc_cost * (mu * conc_values)
                + self._proc_cost * abs_values,
                self._weight * conc_values)

class UndeterminedValues(object):
    def __init__(self, processes, species, functions):
        targets = [p.targets for p in processes]
        # read fluxes used to maintain concentration
        # (only if no value was specified!)
        target_species = [s for t in targets for s in t.concentrations
                          if s.value is None]
        indices = [species.ids.index(t.species) for t in target_species];
        self.names = [t.species for t in target_species]
        lb = [s.lower_bound for s in target_species]
        ub = [s.upper_bound for s in target_species]
        self._conc_composition = species.production[:, indices].tocsr()
        self._conc_proc_cost = species.prod_proc_cost[:, indices].tocsr()
        self._weight = species.weight[:, indices].tocsr()
        # read production and degradation fluxes
        # (only if no value was specified!)
        prod_targets = [s for t in targets for s in t.production_fluxes
                        if s.value is None]
        prod_indices = [species.ids.index(t.species) for t in prod_targets];
        deg_targets = [s for t in targets for s in t.degradation_fluxes
                       if s.value is None]
        deg_indices = [species.ids.index(t.species) for t in deg_targets];
        self.names += [t.species for t in chain(prod_targets, deg_targets)]
        lb += [s.lower_bound for s in chain(prod_targets,deg_targets)]
        ub += [s.upper_bound for s in chain(prod_targets,deg_targets)]
        self._composition = hstack([species.production[:,prod_indices],
                                    species.degradation[:,deg_indices]]).tocsr()
        self._proc_cost = hstack([species.prod_proc_cost[:,prod_indices],
                                  species.deg_proc_cost[:,deg_indices]]).tocsr()
        # lower bounds and upper bounds
        self._lb = TargetVector(lb, functions, 0)
        self._ub = TargetVector(ub, functions, 1e5)

    def matrices(self, mu):
        return (hstack([self._conc_composition * mu, self._composition]),
                hstack([self._conc_proc_cost * mu, self._proc_cost]),
                self._weight)
    
    def lb(self, mu):
        return self._lb.compute(mu)

    def ub(self, mu):
        return self._ub.compute(mu)
