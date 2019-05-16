"""Module defining Solver class."""

# python 2/3 compatibility
from __future__ import division, print_function, absolute_import

# global imports
import numpy
from scipy.sparse import coo_matrix, diags, hstack, vstack
import cplex



def is_feasible(lp):
    return lp.solution.get_status() == lp.solution.status.optimal


def is_infeasible(lp):
    flag = lp.solution.get_status()
    return (flag == lp.solution.status.infeasible
            or flag == lp.solution.status.optimal_infeasible)


def unknown_flag_msg(mu, lp):
    return ('At mu = {}: Unknown exit flag {} corresponding to status {}. '
            'Interrupting computation.'.format(
                mu, lp.solution.get_status(), lp.solution.get_status_string()
                ))


class Solver(object):
    """RBA solver."""

    def __init__(self, matrix):
        self.matrix = matrix

    def _store_results(self, mu, lp):
        self.mu_opt = mu
        self.X = numpy.array(lp.solution.get_values())
        self.lambda_ = lp.solution.get_dual_values()
        self._sol_basis = lp.solution.basis.get_basis()

    def solve(self):
        """Compute configuration corresponding to maximal growth rate."""
        self._sol_basis = None
        self.X = self.lambda_ = self.mu_opt = None

        # check that mu=0 is solution
        self.matrix.build_matrices(0)
        lp = self.build_lp()
        lp.solve()
        if is_feasible(lp):
            self._store_results(0.0, lp)
            print('Mu = 0 is feasible.')
        elif is_infeasible(lp):
            print('Mu = 0 is infeasible, check matrix consistency.')
            return
        else:
            print(unknown_flag_msg(0, lp))
            return


        # bissection
        mu_min = 0
        mu_max = 2.5
        mu_test = mu_max
        self._sol_basis = None
        while (mu_max - mu_min) > 1e-6:
            self.matrix.build_matrices(mu_test)
            lp = self.build_lp()
            lp.solve()
            if is_feasible(lp):
                mu_min = mu_test
                self._store_results(mu_test, lp)
                print('Mu = ' + str(mu_test) + ' is feasible.')
            elif is_infeasible(lp):
                mu_max = mu_test
                print('Mu = ' + str(mu_test) + ' is infeasible.')
            else:
                print(unknown_flag_msg(mu_test, lp))
                return
            # next mu to be tested
            mu_test = (mu_min + mu_max) / 2
            
        # Recompute the matrices so that it necessarily correspond to the solution
        self.matrix.build_matrices(self.mu_opt)
        
    def solve_grid(self):
        """Compute configuration corresponding to maximal growth rate."""
        self._sol_basis = None
        self.X = self.lambda_ = self.mu_opt = None

        # check that mu=0 is solution
        self.matrix.build_matrices(0)
        lp = self.build_lp()
        lp.solve()
        if is_feasible(lp):
            self._store_results(0.0, lp)
        elif is_infeasible(lp):
            print('Mu = 0 is infeasible, check matrix consistency.')
            return
        else:
            print(unknown_flag_msg(0, lp))
            return

        # grid
        vec_mu = numpy.arange(0, 0.8, 0.001)
        for mu_test in vec_mu:
            self.matrix.build_matrices(mu_test)
            lp = self.build_lp()
            lp.solve()
            if is_feasible(lp):
                print('Mu = ' + str(mu_test) + ' is feasible.')
                self._store_results(mu_test, lp)
            elif is_infeasible(lp):
                print('Mu = ' + str(mu_test) + ' is infeasible.')
            else:
                print(unknown_flag_msg(mu_test, lp))
                return

    def build_lp(self):
        """
        Build CPLEX problem based on current matrices.

        Returns
        -------
        cplex object
            LP problem with default parameters corresponding to current
            matrices.

        """
        # preprocess matrices
        # rescale concentration columns?
        lhs = self.matrix.A
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
        lp_problem.variables.add(names=self.matrix.col_names)
        lp_problem.variables.set_lower_bounds(zip(self.matrix.col_names, self.matrix.LB))
        lp_problem.variables.set_upper_bounds(zip(self.matrix.col_names, self.matrix.UB))
        lp_problem.objective.set_sense(lp_problem.objective.sense.minimize)
        lp_problem.objective.set_linear(zip(self.matrix.col_names, self.matrix.f))

        lp_problem.linear_constraints.add(names=self.matrix.row_names)
        lp_problem.linear_constraints.set_linear_components(zip(self.matrix.row_names, rows))
        lp_problem.linear_constraints.set_rhs(zip(self.matrix.row_names, self.matrix.b))
        lp_problem.linear_constraints.set_senses(zip(self.matrix.row_names,self.matrix.row_signs))
        # set starting point (not exactly sure how this works)
        if self._sol_basis is not None:
            lp_problem.start.set_start(
                self._sol_basis[0], self._sol_basis[1], self.X, [], [],
                self.lambda_)
        return lp_problem
