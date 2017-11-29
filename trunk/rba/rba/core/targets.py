"""Module processing target information."""

# python 2/3 compatibility
from __future__ import division, print_function, absolute_import

# global imports
import numpy
from collections import namedtuple
from scipy.sparse import hstack

# local imports
from rba.core.parameter_vector import ParameterVector

# class used to communicate target information
TargetSet = namedtuple('TargetSet', 'names values lb ub composition '
                       'processing_cost weight')


class Targets(object):
    """
    Class computing target-related substructures.

    Attributes
    ----------
    determined_targets : rba.core.targets.DeterminedTargets
        Target fluxes for metabolites (only fixed values).
    determined_targets : rba.core.targets.UndeterminedTargets
        Undetermined flux values (only bounds were given).
    target_reactions : rba.core.targets.TargetReactions

    """

    def __init__(self, targets, species, parameters):
        """
        Constructor.

        Parameters
        ----------
        targets : rba.xml.RbaTargets
            Structure containing target information.
        species : rba.core.species.Species
            Species information.
        parameters : rba.core.parameters.Parameters
            Parameter information

        """
        target_groups = targets.target_groups
        # extract target values (absolute + concentration related fluxes)
        self.determined_targets = DeterminedTargets(
            target_groups, species, parameters
            )
        self.undetermined_targets = UndeterminedTargets(
            target_groups, species, parameters
            )
        # extract target reactions
        self.target_reactions = TargetReactions(target_groups, parameters)


class DeterminedTargets(object):
    """Class computing metabolite target production/consumption."""

    def __init__(self, targets, species, parameters):
        """
        Constructor.

        Parameters
        ----------
        targets : rba.xml.ListOfTargetGroups
            xml structure containing target information.
        species : rba.core.species.Species
            Species information.
        parameters : rba.core.parameters.Parameters
            Parameter information.

        """
        # extract target substructures
        conc_targets = []
        prod_targets = []
        deg_targets = []
        for target_group in targets:
            conc_targets += [t for t in target_group.concentrations
                             if t.value is not None]
            prod_targets += [t for t in target_group.production_fluxes
                             if t.value is not None]
            deg_targets += [t for t in target_group.degradation_fluxes
                            if t.value is not None]
        # extract fluxes used to maintain concentration
        conc_set = extract_targets(conc_targets, species)
        self._conc_composition = conc_set.composition
        self._conc_proc_cost = conc_set.processing_cost
        self._weight = conc_set.weight
        self._conc_values = ParameterVector(conc_set.values, parameters)
        # extract absolute production and degradation fluxes
        prod_set = extract_targets(prod_targets, species)
        deg_set = extract_targets(deg_targets, species, degradation=True)
        self._values = ParameterVector(prod_set.values + deg_set.values,
                                       parameters)
        self._composition = hstack([prod_set.composition,
                                    deg_set.composition]).tocsr()
        self._proc_cost = hstack([prod_set.processing_cost,
                                  deg_set.processing_cost]).tocsr()

    def compute(self, mu):
        """Compute composition, processing cost and weight of targets."""
        conc_values = self._conc_values.compute()
        abs_values = self._values.compute()
        return (self._conc_composition * (mu * conc_values)
                + self._composition * abs_values,
                self._conc_proc_cost * (mu * conc_values)
                + self._proc_cost * abs_values,
                self._weight * conc_values)


class UndeterminedTargets(object):
    """
    Class computing upper/lower bounds of process-related fluxes.

    Attributes
    ----------
    names : list of str
        Names of metabolites with undetermined flux.

    """

    def __init__(self, targets, species, parameters):
        """
        Constructor.

        Parameters
        ----------
        targets : rba.xml.ListOfTargetGroups
            xml structure containing target information.
        species : rba.core.species.Species
            Species information.
        parameters : rba.core.parameters.Parameters
            Parameter information.

        """
        # extract target substructures
        conc_targets = []
        prod_targets = []
        deg_targets = []
        for target_group in targets:
            conc_targets += [t for t in target_group.concentrations
                             if t.value is None]
            prod_targets += [t for t in target_group.production_fluxes
                             if t.value is None]
            deg_targets += [t for t in target_group.degradation_fluxes
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
        self._lb = ParameterVector(conc_set.lb + prod_set.lb + deg_set.lb,
                                   parameters)
        self._ub = ParameterVector(conc_set.ub + prod_set.ub + deg_set.ub,
                                   parameters)
        self.f = numpy.zeros(len(conc_targets) + len(prod_targets)
                             + len(deg_targets))

    def matrices(self, mu):
        """
        Return composition, processing cost and weight matrices of targets.

        Parameters
        ----------
        mu : float
            growth rate.

        Returns
        -------
        Tuple of matrices : First element is the composition matrix, second
            element is the processing cost matrix and last element is the
            weight matrix.

        """
        return (hstack([self._conc_composition * mu, self._composition]),
                hstack([self._conc_proc_cost * mu, self._proc_cost]),
                self._weight)

    def lb(self):
        """
        Return lower bound vector.

        Returns
        -------
        Vector containing lower bounds on target fluxes.

        """
        return self._lb.compute()

    def ub(self):
        """
        Return upper bound vector.

        Returns
        -------
        Vector containing upper bounds on target fluxes.

        """
        return self._ub.compute()


class TargetReactions(object):
    """
    Class computing process-related reaction fluxes.

    Attributes
    ----------
        value_reactions: list of reactions whose flux is determined by a
            process.
        lb_reactions: list of reactions whose lower bound is determined by
            a process.
        ub_reactions: list of reactions whose upper bound is determined by
            a process.

    """

    def __init__(self, targets, parameters):
        """
        Constructor.

        Parameters
        ----------
        targets : rba.xml.ListOfTargetGroups
            xml structure containing target information.
        species : rba.core.species.Species
            Species information.
        parameters : rba.core.parameters.Parameters
            Parameters.

        """
        self.value_reactions = []
        self.lb_reactions = []
        self.ub_reactions = []
        val = []
        lb = []
        ub = []
        for target_group in targets:
            for target in target_group.reaction_fluxes:
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
        self._value = ParameterVector(val, parameters)
        self._lb = ParameterVector(lb, parameters)
        self._ub = ParameterVector(ub, parameters)

    def value(self):
        """
        Return vector of reaction fluxes.

        Returns
        -------
        Vector containing flux values for reactions whose flux is
        determined by a process.

        """
        return self._value.compute()

    def lb(self):
        """
        Return vector of lower bounds.

        Returns
        -------
        Vector containing lower boundvalues for reactions whose
        lower bound is determined by a process.

        """
        return self._lb.compute()

    def ub(self):
        """
        Return vector of upper bounds.

        Returns
        -------
        Vector containing upper bound values for reactions whose
        upper bound is determined by a process.

        """
        return self._ub.compute()


def extract_targets(targets, species, degradation=False):
    """
    Extract basic target information for fixed targets.

    Parameters
    ----------
    targets : list of rba.xml.TargetSpecies
        Target species to extract.
    species : rba.core.species.Species
        Species information.
    degradation : bool
        If set to False, composition and processing cost
        are computed assuming that targets are being produced. If set to
        True, composition and processing cost are computed assuming targets
        are being degraded.

    Returns
    -------
    TargetSet : object containing names of targets, target values,
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
