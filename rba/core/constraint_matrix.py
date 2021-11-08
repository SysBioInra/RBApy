"""Module defining ConstraintMatrix class."""

# python 2/3 compatibility
from __future__ import division, print_function, absolute_import

# global imports
import numpy
from scipy.sparse import coo_matrix, diags, hstack, vstack

# local imports
from rba.core.constraint_blocks import ConstraintBlocks


class ConstraintMatrix(object):
    """
    Class building constraint matrix.

    Attributes:
        col_names: Linear problem column names (decision variables).
        reaction_cols: Indices of columns corresponding to reactions.
        enzyme_cols: Indices of columns corresponding to enzymes.
        process_cols: Indices of columns corresponding to processes.
        target_cols: Indices of columns corresponding to targets.
        row_names: Linear problem row names (constraints).
        row_signs: Linear problem row signs (equality or inequality).
        UB: Linear problem upper bounds.
        LB: Linear problem lower bounds.
        f: Linear problem objective function.
        A: Linear problem matrix (left-hand side).
        b: Linear problem right-hand side.

    """

    def __init__(self, model):
        """
        Build constraint matrix from model.

        Parameters
        ----------
        model : rba.RbaModel
            RBA model.

        """
        self._blocks = ConstraintBlocks(model)
        # convenience variables
        reactions = self._blocks.metabolism.reactions
        enzymes = self._blocks.enzymes.ids
        processes = self._blocks.processes.ids
        undetermined_fluxes = self._blocks.targets.undetermined_targets.names
        compartments = self._blocks.density.compartments
        nb_reactions = len(reactions)
        nb_enzymes = len(enzymes)
        nb_processes = len(processes)
        nb_undetermined = len(undetermined_fluxes)
        nb_compartments = len(compartments)
        # column information
        self.col_names = (reactions
                          + [e for e in enzymes]
                          + [p + '_machinery' for p in processes]
                          + [m + '_target_flux' for m in undetermined_fluxes])
        self.reaction_cols = numpy.arange(nb_reactions)
        self.enzyme_cols = nb_reactions + numpy.arange(nb_enzymes)
        self.process_cols = (nb_reactions + nb_enzymes +
                             numpy.arange(nb_processes))
        self.target_cols = (nb_reactions + nb_enzymes + nb_processes +
                            numpy.arange(nb_undetermined))
        # row information
        self.row_names = (self._blocks.metabolism.internal
                          + [p + '_capacity' for p in processes]
                          + [e + '_forward_capacity' for e in enzymes]
                          + [e + '_backward_capacity' for e in enzymes]
                          + [c + '_density' for c in compartments])
        self.row_signs = (['E'] * len(self._blocks.metabolism.internal)
                          + self._blocks.processes.capacity_signs
                          + ['L'] * 2 * nb_enzymes
                          + self._blocks.density.signs)
        # constant building blocks
        self._empty_ExPU = coo_matrix((nb_enzymes,
                                       nb_processes + nb_undetermined))
        self._empty_PxR = coo_matrix((nb_processes, nb_reactions))
        self._empty_CxR = coo_matrix((nb_compartments, nb_reactions))
        self._empty_2E = numpy.zeros(2 * nb_enzymes)
        # indicator matrices
        R_ind = [reactions.index(r)
                 for r in self._blocks.enzymes.reaction_catalyzed]
        self._R_to_E = coo_matrix(([1]*nb_enzymes, (range(nb_enzymes), R_ind)),
                                  shape=(nb_enzymes, nb_reactions))
        target_reactions = self._blocks.targets.target_reactions
        self._value_reaction_cols = self.reaction_cols[
            [reactions.index(r) for r in target_reactions.value_reactions]
            ]
        self._lb_reaction_cols = self.reaction_cols[
            [reactions.index(r) for r in target_reactions.lb_reactions]
            ]
        self._ub_reaction_cols = self.reaction_cols[
            [reactions.index(r) for r in target_reactions.ub_reactions]
            ]
        # set remaining attributes to None
        self.A = self.b = self.LB = self.UB = self.f = None

    def build_matrices(self, mu):
        """
        Build LP matrices corresponding to given growth-rate.

        Args:
            mu: growth_rate
        """
        # update parameters
        self._blocks.parameters.update_growth_rate(mu)

        # build A
        enzymes = self._blocks.enzymes
        processes = self._blocks.processes
        targets = self._blocks.targets
        density = self._blocks.density
        # mu-dependent blocks
        u_composition, u_proc_cost, u_weight \
            = targets.undetermined_targets.matrices(mu)
        process_capacity = processes.capacity.compute()
        forward, backward = enzymes.efficiencies()
        # stoichiometry constraints
        metab_rows = hstack([self._blocks.metabolism.S,
                             mu * enzymes.machinery.composition,
                             mu * processes.machinery.composition,
                             u_composition])
        # capacity constraints
        process_rows = hstack([self._empty_PxR,
                               mu * enzymes.machinery.processing_cost,
                               mu * processes.machinery.processing_cost
                               - diags(process_capacity),
                               u_proc_cost])
        forward_rows = hstack(
            [self._R_to_E, -diags(forward), self._empty_ExPU]
            )
        backward_rows = hstack(
            [-self._R_to_E, -diags(backward), self._empty_ExPU]
            )
        # density constraints
        c_indices = density.compartment_indices
        density_rows = hstack([self._empty_CxR,
                               enzymes.machinery.weight[c_indices],
                               processes.machinery.weight[c_indices],
                               u_weight[c_indices]])
        self.A = vstack([metab_rows, process_rows,
                         forward_rows, backward_rows, density_rows])

        # build b
        # gather mu-dependent blocks
        fluxes, processing, weight = targets.determined_targets.compute(mu)
        density_rows = density.values.compute() - weight[c_indices].T
        # build vector
        self.b = numpy.concatenate([-fluxes, -processing,
                                    self._empty_2E, density_rows])

        # update lower bounds and upper bounds
        self.LB = numpy.concatenate([self._blocks.metabolism.lb(),
                                     self._blocks.enzymes.lb,
                                     processes.lb,
                                     targets.undetermined_targets.lb()])
        self.UB = numpy.concatenate([self._blocks.metabolism.ub(),
                                     self._blocks.enzymes.ub,
                                     processes.ub,
                                     targets.undetermined_targets.ub()])
        self.f = numpy.concatenate([self._blocks.metabolism.f,
                                    self._blocks.enzymes.f,
                                    processes.f,
                                    targets.undetermined_targets.f])
        # target reactions
        self.LB[self._lb_reaction_cols] = targets.target_reactions.lb()
        self.UB[self._ub_reaction_cols] = targets.target_reactions.ub()
        r_fluxes = targets.target_reactions.value()
        self.LB[self._value_reaction_cols] = r_fluxes
        self.UB[self._value_reaction_cols] = r_fluxes

    def set_medium(self, medium):
        """
        Change medium composition.

        Args:
            medium: dict mapping metabolite prefixes with their concentration.
        """
        self._blocks.set_medium(medium)
