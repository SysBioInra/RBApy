"""Module defining RbaSolver class."""

# python 2/3 compatibility
from __future__ import division, print_function

# global imports
import numpy
from scipy.sparse import coo_matrix, diags, hstack, vstack
import cplex


class RbaSolver(object):
    """
    Class solving RBA problem.

    Attributes:
        col_names: Linear problem column names (decision variables).
        reaction_cols: Indices of columns corresponding to reactions.
        enzyme_cols: Indices of columns corresponding to enzymes.
        process_cols: Indices of columns corresponding to processes.
        species_cols: Indices of columns corresponding to enzymes.
        row_names: Linear problem row names (constraints).
        row_signs: Linear problem row signs (equality or inequality).
        UB: Linear problem upper bounds.
        LB: Linear problem lower bounds.
        f: Linear problem objective function.
        A: Linear problem matrix (left-hand side).
        b: Linear problem right-hand side.
        X: Solution to RBA problem (if one was found).
        lambda_: Dual solution to RBA problem (if one was found).
        mu_opt: Growth rate at solution (if one was found).

    """

    def __init__(self, blocks):
        """
        Constructor.

        Args:
            blocks (RbaMatrices): RBA substructures used to compute matrices.
        """
        self._blocks = blocks
        # convenience variables
        reactions = blocks.metabolism.reactions
        enzymes = blocks.enzymes.ids
        processes = blocks.processes.ids
        undetermined_fluxes = blocks.processes.undetermined_values.names
        compartments = blocks.density.compartments
        nb_reactions = len(reactions)
        nb_enzymes = len(enzymes)
        nb_processes = len(processes)
        nb_undetermined = len(undetermined_fluxes)
        nb_compartments = len(compartments)
        # column information
        self.col_names = (reactions + enzymes
                          + [p + '_machinery' for p in processes]
                          + [m + '_flux' for m in undetermined_fluxes])
        self.reaction_cols = numpy.arange(nb_reactions)
        self.enzyme_cols = (self.reaction_cols[-1] + 1
                            + numpy.arange(nb_enzymes))
        self.process_cols = (self.enzyme_cols[-1] + 1
                             + numpy.arange(nb_processes))
        self.species_cols = (self.process_cols[-1] + 1
                             + numpy.arange(nb_undetermined))
        # row information
        self.row_names = (blocks.metabolism.internal
                          + [p + '_capacity' for p in processes]
                          + [e + '_forward_capacity' for e in enzymes]
                          + [e + '_backward_capacity' for e in enzymes]
                          + [c + '_density' for c in compartments])
        self.row_signs = (['E'] * len(blocks.metabolism.internal)
                          + blocks.processes.capacity_signs
                          + ['L'] * 2 * nb_enzymes
                          + blocks.density.signs)
        # upper bound, lower bound, objective function
        nb_cols = len(self.col_names)
        self.UB = 1e5 * numpy.ones(nb_cols)
        self.LB = numpy.zeros(nb_cols)
        self.f = numpy.zeros(nb_cols)
        self.LB[self.reaction_cols] = (
            -1e3 * numpy.array(blocks.metabolism.reversibility)
            )
        self.LB[self.reaction_cols[blocks.metabolism.zero_lb]] = 0
        self.UB[self.reaction_cols[blocks.metabolism.zero_ub]] = 0
        self.f[self.enzyme_cols] = 1
        # constant building blocks
        self._empty_ExPU = coo_matrix((nb_enzymes,
                                       nb_processes + nb_undetermined))
        self._empty_PxR = coo_matrix((nb_processes, nb_reactions))
        self._empty_CxR = coo_matrix((nb_compartments, nb_reactions))
        self._empty_2E = numpy.zeros(2 * nb_enzymes)
        # indicator matrices
        R_ind = [reactions.index(r) for r in blocks.enzymes.reaction_catalyzed]
        self._R_to_E = coo_matrix(([1]*nb_enzymes, (range(nb_enzymes), R_ind)),
                                  shape=(nb_enzymes, nb_reactions))
        target_reactions = blocks.processes.target_reactions
        self._value_reaction_cols \
            = self.reaction_cols[[reactions.index(r)
                                  for r in target_reactions.value_reactions]]
        self._lb_reaction_cols \
            = self.reaction_cols[[reactions.index(r)
                                  for r in target_reactions.lb_reactions]]
        self._ub_reaction_cols \
            = self.reaction_cols[[reactions.index(r)
                                  for r in target_reactions.ub_reactions]]
        # set remaining attributes to None
        self.A = self.b = None
        self.mu_opt = self.X = self.lambda_ = self._sol_basis = None

    def build_matrices(self, mu):
        """
        Build LP matrices corresponding to given growth-rate.

        Args:
            mu: growth_rate
        """
        # build A
        enzymes = self._blocks.enzymes
        processes = self._blocks.processes
        density = self._blocks.density
        # mu-dependent data
        u_composition, u_proc_cost, u_weight \
            = self._blocks.processes.undetermined_values.matrices(mu)
        process_capacity = processes.capacity.compute(mu)
        (forward, backward) = enzymes.efficiency.compute(mu)
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
        # gather mu-dependent data
        (fluxes, processing, weight) = processes.target_values.compute(mu)
        density_rows = density.values.compute(mu) - weight[c_indices].T
        # build vector
        self.b = numpy.concatenate([-fluxes, -processing,
                                    self._empty_2E, density_rows])

        # update lower bounds and upper bounds
        # undetermined metabolites
        self.LB[self.species_cols] = processes.undetermined_values.lb(mu)
        self.UB[self.species_cols] = processes.undetermined_values.ub(mu)
        # target reactions
        self.LB[self._lb_reaction_cols] = processes.target_reactions.lb(mu)
        self.UB[self._ub_reaction_cols] = processes.target_reactions.ub(mu)
        r_fluxes = self._blocks.processes.target_reactions.value(mu)
        self.LB[self._value_reaction_cols] = r_fluxes
        self.UB[self._value_reaction_cols] = r_fluxes

    def solve(self):
        """Compute configuration corresponding to maximal growth rate."""
        self._sol_basis = None

        # check that mu=0 is solution
        self.build_matrices(0)
        lp = self.build_lp()
        lp.solve()
        exit_flag = lp.solution.get_status()
        if exit_flag == lp.solution.status.optimal:
            self.X = numpy.array(lp.solution.get_values())
            self.lambda_ = lp.solution.get_dual_values()
            self.mu_opt = 0.0
            self._sol_basis = lp.solution.basis.get_basis()
        elif (exit_flag == lp.solution.status.infeasible
              or (exit_flag == lp.solution.status.optimal_infeasible)):
            print('Mu = 0 is infeasible, check matrix consistency')
            return
        else:
            print('At mu = 0: Unknown exit flag {} corresponding to '
                  'status {}. Interrupting computation.'
                  .format(exit_flag, lp.solution.get_status_string()))
            return

        # bissection
        mu_min = 0.0
        mu_max = 2.5
        mu_test = mu_max
        while (mu_max - mu_min) > 1e-4:
            self.build_matrices(mu_test)
            lp = self.build_lp()
            lp.solve()
            exit_flag = lp.solution.get_status()
            print(mu_test, lp.solution.get_status_string())
            if exit_flag == lp.solution.status.optimal:
                mu_min = mu_test
                self.X = numpy.array(lp.solution.get_values())
                self.lambda_ = lp.solution.get_dual_values()
                self.mu_opt = mu_min
                self._sol_basis = lp.solution.basis.get_basis()
            elif (exit_flag == lp.solution.status.infeasible
                  or exit_flag == lp.solution.status.optimal_infeasible):
                mu_max = mu_test
            else:
                print('At mu = {}: Unknown exit flag {} corresponding to '
                      'status {}. Interrupting computation.'
                      .format(mu_test, exit_flag,
                              lp.solution.get_status_string()))
                return
            # next mu to be tested
            mu_test = (mu_min + mu_max) / 2
        print(mu_min)

    def build_lp(self):
        """
        Build CPLEX problem based on current matrices.

        Returns:
            cplex object with default parameters corresponding to current
            matrices.

        """
        # preprocess matrices
        # rescale concentration columns?
        lhs = self.A
        # scaling_factor = 1000
        # scaling = numpy.ones(lhs.shape[1])
        # scaling[numpy.concatenate([self.enzyme_cols, self.process_cols,
        # self.species_cols])] \
        # = 1.0/scaling_factor
        # lhs *= diags(scaling)

        # transform inequality and equality constraints to CPLEX row format
        lhs = lhs.tolil()
        rows = []
        for nz_ind, data in zip(lhs.rows, lhs.data):
            rows.append(cplex.SparsePair(nz_ind, data))

        # define problem
        lp_problem = cplex.Cplex()
        # set parameters
        lp_problem.objective.set_sense(lp_problem.objective.sense.minimize)
        lp_problem.parameters.feasopt.tolerance.set(1e-9)
        lp_problem.parameters.simplex.tolerances.feasibility.set(1e-9)
        lp_problem.parameters.simplex.tolerances.optimality.set(1e-9)
        lp_problem.parameters.simplex.tolerances.markowitz.set(0.1)
        lp_problem.parameters.barrier.convergetol.set(1e-9)
        # agressive scaling
        lp_problem.parameters.read.scale.set(1)
        # Threads: the default (0) means that Cplex decides automatically
        # how many threads to use
        # lp_problem.parameters.threads.set(0)
        lp_problem.set_results_stream(None)
        # define columns and add rows
        lp_problem.variables.add(
            obj=self.f, ub=self.UB, lb=self.LB, names=self.col_names
            )
        lp_problem.linear_constraints.add(
            lin_expr=rows, rhs=self.b, senses=self.row_signs,
            names=self.row_names
            )
        # set starting point (not exactly sure how this works)
        if self._sol_basis is not None:
            lp_problem.start.set_basis(self._sol_basis[0], self._sol_basis[1])
            # lp_problem.start.set_start(
            # self._sol_basis[0], self._sol_basis[1], self.X, [], [],
            # self.lambda_)
        return lp_problem
