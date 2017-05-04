
from scipy.sparse import hstack, csr_matrix, lil_matrix

from target_vector import TargetVector

class Processes(object):
    def __init__(self, processes, species, functions):
        ## extract ids
        self.ids = [p.id for p in processes]

        ## extract machinery related information
        # extract machinery
        machinery = [p.machinery.machinery_composition for p in processes]
        self.machinery = species.create_machinery(machinery)
        # extract capacities
        self.capacity = TargetVector([p.machinery.capacity for p in processes],
                                     functions)
        
        ## extract targets
        # extract target values (absolute + concentration related fluxes)
        self.target_values = TargetValues(processes, species, functions)
        # extract target reactions
        target_reactions = [r for p in processes for r in p.targets.reaction_fluxes]
        self.target_reaction = [t.reaction for t in target_reactions]
        self.reaction_value = TargetVector(target_reactions, functions)
        # extract lower bounds
        target_reactions = [r for p in processes for r in p.targets.lower_bounds]
        self.lb_reaction = [t.reaction for t in target_reactions]
        self.lb = TargetVector(target_reactions, functions)
        # extract upper bounds
        target_reactions = [r for p in processes for r in p.targets.upper_bounds]
        self.ub_reaction = [t.reaction for t in target_reactions]
        self.ub = TargetVector(target_reactions, functions)

class TargetValues(object):
    def __init__(self, processes, species, functions):
        targets = [p.targets for p in processes]
        # read fluxes used to maintain concentration
        target_species = [s for t in targets for s in t.concentrations]
        indices = [species.ids.index(t.species) for t in target_species];
        self._conc_values = TargetVector(target_species, functions)
        self._conc_composition = species.production[:, indices].tocsr()
        self._conc_proc_cost = species.prod_proc_cost[:, indices].tocsr()
        self._weight = species.weight[:, indices].tocsr()
        # read production and degradation fluxes
        prod_targets = [s for t in targets for s in t.production_fluxes]
        prod_indices = [species.ids.index(t.species) for t in prod_targets];
        deg_targets = [s for t in targets for s in t.degradation_fluxes]
        deg_indices = [species.ids.index(t.species) for t in deg_targets];
        self._values = TargetVector(prod_targets+deg_targets, functions)
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
