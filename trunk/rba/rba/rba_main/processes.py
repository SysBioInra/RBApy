"""
Module defining Processes, TargetValues, UndeterminedValues, 
TargetReactions classes.
"""

# python 2/3 compatibility
from __future__ import division, print_function

# global imports
from collections import namedtuple
from scipy.sparse import hstack

# local imports
from rba.rba_main.target_vector import TargetVector

# class used to communicate target information
TargetSet = namedtuple('TargetSet', 'names values lb ub composition '
                       'processing_cost weight')

class Processes(object):
    """
    Class computing process-related substructures.

    Attributes:
        ids: list of process identifiers.
        machinery (list of species:Machinery): list of objects containing all
            composition-related information of process machineries.
        capacity (target_vector.TargetVector): object containing right-hand
            side of process capacity constraints.
        capacity_signs: list of signs of capacity constraints ('E' for equality,
            'L' for inequality).
        target_values (processes.TargetValues): object used to compute target
            fluxes for metabolites (only fixed values).
        undetermined_values (processes.UndeterminedValues): object storing
            undetermined flux values (only lower bound and/or upper bound
            was given).
    """

    def __init__(self, processes, species, functions):
        """
        Constructor.

        Args:
            processes: xml structure containing process information.
            species (species.Species): object containing information about
                species in the system.
            functions: dict mapping known function ids with objects computing
                them.
        """
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
        for process in processes:
            if process.machinery.capacity.value is not None:
                values.append(process.machinery.capacity.value)
                self.capacity_signs.append('E')
            elif process.machinery.capacity.upper_bound is not None:
                values.append(process.machinery.capacity.value)
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
        self.target_reactions = TargetReactions(processes, functions)


class TargetValues(object):
    """
    Class computing metabolite fluxes produced/consumed by processes.
    """

    def __init__(self, processes, species, functions):
        """
        Constructor.

        Args:
            processes: xml structure containing process information.
            species (species.Species): object containing information about
                species in the system.
            functions: dict mapping known function ids with objects computing
                them.
        """
        targets = [p.targets for p in processes]
        # extract target substructures
        conc_targets = []
        prod_targets = []
        deg_targets = []
        for target in targets:
            conc_targets += [t for t in target.concentrations
                             if t.value is not None]
            prod_targets += [t for t in target.production_fluxes
                             if t.value is not None]
            deg_targets += [t for t in target.degradation_fluxes
                            if t.value is not None]
        # extract fluxes used to maintain concentration
        conc_set = extract_targets(conc_targets, species)
        self._conc_composition = conc_set.composition
        self._conc_proc_cost = conc_set.processing_cost
        self._weight = conc_set.weight
        self._conc_values = TargetVector(conc_set.values, functions)
        # extract absolute production and degradation fluxes
        prod_set = extract_targets(prod_targets, species)
        deg_set = extract_targets(deg_targets, species, degradation=True)
        self._values = TargetVector(prod_set.values + deg_set.values, functions)
        self._composition = hstack([prod_set.composition,
                                    deg_set.composition]).tocsr()
        self._proc_cost = hstack([prod_set.processing_cost,
                                  deg_set.processing_cost]).tocsr()

    def compute(self, mu):
        """
        Compute composition, processing cost and weight of targets.

        Args:
            mu: growth rate.
        """
        conc_values = self._conc_values.compute(mu)
        abs_values = self._values.compute(mu)
        return (self._conc_composition * (mu * conc_values)
                + self._composition * abs_values,
                self._conc_proc_cost * (mu * conc_values)
                + self._proc_cost * abs_values,
                self._weight * conc_values)


class UndeterminedValues(object):
    """
    Class computing upper/lower bounds of process-related fluxes.

    Attributes:
        names: names of metabolites which flux is undetermined.
    """

    def __init__(self, processes, species, functions):
        """
        Constructor.

        Args:
            processes: xml structure containing process information.
            species (species.Species): object containing information about
                species in the system.
            functions: dict mapping known function ids with objects computing
                them.
        """
        targets = [p.targets for p in processes]
        # extract target substructures
        conc_targets = []
        prod_targets = []
        deg_targets = []
        for target in targets:
            conc_targets += [t for t in target.concentrations
                             if t.value is None]
            prod_targets += [t for t in target.production_fluxes
                             if t.value is None]
            deg_targets += [t for t in target.degradation_fluxes
                            if t.value is None]
        # extract fluxes used to maintain concentration
        conc_set = extract_targets(conc_targets, species)
        # read production and degradation fluxes
        prod_set = extract_targets(prod_targets, species)
        deg_set = extract_targets(deg_targets, species, degradation=True)
        self.names = conc_set.names + prod_set.names + deg_set.names
        self._conc_composition = conc_set.composition
        self._conc_proc_cost = conc_set.processing_cost
        self._weight = conc_set.weight
        self._composition = hstack([prod_set.composition,
                                    deg_set.composition]).tocsr()
        self._proc_cost = hstack([prod_set.processing_cost,
                                  deg_set.processing_cost]).tocsr()
        self._lb = TargetVector(conc_set.lb + prod_set.lb + deg_set.lb,
                                functions, 0)
        self._ub = TargetVector(conc_set.ub + prod_set.ub + deg_set.ub,
                                functions, 1e5)

    def matrices(self, mu):
        """
        Return composition, processing cost and weight matrices of targets.

        Args:
            mu: growth rate.

        Returns:
            Tuple of matrices. First element is the composition matrix, second
            element is the processing cost matrix and last element is the
            weight matrix.
        """
        return (hstack([self._conc_composition * mu, self._composition]),
                hstack([self._conc_proc_cost * mu, self._proc_cost]),
                self._weight)

    def lb(self, mu):
        """
        Return lower bound vector.

        Args:
            mu: growth rate.

        Returns:
            Vector containing lower bounds on target fluxes.
        """
        return self._lb.compute(mu)

    def ub(self, mu):
        """
        Return upper bound vector.

        Args:
            mu: growth rate.

        Returns:
            Vector containing upper bounds on target fluxes.
        """
        return self._ub.compute(mu)


class TargetReactions(object):
    """
    Class computing process-related reaction fluxes.

    Attributes:
        value_reactions: list of reactions whose flux is determined by a
            process.
        lb_reactions: list of reactions whose lower bound is determined by
            a process.
        ub_reactions: list of reactions whose upper bound is determined by
            a process.
    """

    def __init__(self, processes, functions):
        """
        Constructor.

        Args:
            processes: xml structure containing process information.
            species (species.Species): object containing information about
                species in the system.
            functions: dict mapping known function ids with objects computing
                them.
        """
        self.value_reactions = []
        self.lb_reactions = []
        self.ub_reactions = []
        val = []
        lb = []
        ub = []
        for process in processes:
            for target in process.targets.reaction_fluxes:
                if target.value is not None:
                    val.append(target.value)
                    self.value_reactions.append(target.reaction)
                else:
                    if target.lower_bound is not None:
                        lb.append(target.lower_bound)
                        self.lb_reactions.append(target.reaction)
                    if target.upper_bound is not None:
                        ub.append(target.upper_bound)
                        self.ub_reactions.append(target.reaction)
        self._value = TargetVector(val, functions)
        self._lb = TargetVector(lb, functions)
        self._ub = TargetVector(ub, functions)

    def value(self, mu):
        """
        Return vector of reaction fluxes.

        Args:
            mu: growth rate

        Returns:
            Vector containing flux values for reactions whose flux is determined
            by a process.
        """
        return self._value.compute(mu)

    def lb(self, mu):
        """
        Return vector of lower bounds.

        Args:
            mu: growth rate

        Returns:
            Vector containing lower boundvalues for reactions whose
            lower bound is determined by a process.
        """
        return self._lb.compute(mu)

    def ub(self, mu):
        """
        Return vector of upper bounds.

        Args:
            mu: growth rate

        Returns:
            Vector containing upper bound values for reactions whose
            upper bound is determined by a process.
        """
        return self._ub.compute(mu)

def extract_targets(targets, species, degradation=False):
    """
    Extract basic target information for fixed targets.

    Args:
        targets: list of xml structure containing target species.
        species (species.Species): object containing information about
            species in the system.
        degradation: flag. If set to False, composition and processing cost
            are computed assuming that targets are being produced. If set to
            True, composition and processing cost are computed assuming targets
            are being degraded.

    Returns:
        TargetSet object containing names of targets, target values,
        lower bounds, upper bounds, composition matrix, processing cost matrix
        and weight vector.
    """
    names = [t.species for t in targets]
    indices = [species.ids.index(n) for n in names]
    values = [t.value for t in targets]
    lb = [t.lower_bound for t in targets]
    ub = [t.upper_bound for t in targets]
    weight = species.weight[:, indices].tocsr()
    if not degradation:
        composition = species.production[:, indices].tocsr()
        proc_cost = species.prod_proc_cost[:, indices].tocsr()
    else:
        composition = species.degradation[:, indices].tocsr()
        proc_cost = species.deg_proc_cost[:, indices].tocsr()
    return TargetSet(names, values, lb, ub, composition, proc_cost, weight)
