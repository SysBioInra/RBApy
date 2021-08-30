"""Module defining Solver class."""

# python 2/3 compatibility
from __future__ import division, print_function, absolute_import

# global imports
import abc
import numpy
# from scipy.sparse import diags
import sys
try:
    import cplex
except ModuleNotFoundError:
    cplex = None


class Solver(object):
    """RBA solver."""

    def __init__(self, matrix, lp_solver=None, bissection_tol=1e-6, max_bissection_iters=None, verbose=False):
        self.matrix = matrix

        if (cplex and lp_solver is None) or lp_solver == 'cplex':
            self.lp_solver = CplexLpSolver(self)
        elif lp_solver == 'cplex_optlang':
            self.lp_solver = OptlangLpSolver(self, 'cplex')
        elif (not cplex and lp_solver is None) or lp_solver == 'glpk':
            self.lp_solver = OptlangLpSolver(self, 'glpk')
        elif lp_solver == 'glpk_exact':
            self.lp_solver = OptlangLpSolver(self, 'glpk_exact')
        elif lp_solver == 'scipy':
            self.lp_solver = OptlangLpSolver(self, 'scipy')
        else:
            raise NotImplementedError('LP solver `{}` is not implemented'.format(lp_solver))

        self.bissection_tol = bissection_tol
        self.max_bissection_iters = max_bissection_iters

        self.verbose = verbose

    def solve(self):
        """Compute configuration corresponding to maximal growth rate."""
        self._sol_basis = None
        self.X = self.lambda_ = self.mu_opt = None

        # check that μ=0 is solution
        if self.verbose:
            print('  Checking μ = 0 is feasible ...', end='')
            sys.stdout.flush()
        self.matrix.build_matrices(0)
        self.lp_solver.build_lp()
        self.lp_solver.solve_lp()
        if self.lp_solver.is_feasible():
            self.lp_solver.store_results(0.0)
            if self.verbose:
                print(' μ = 0 is feasible.')
        elif self.lp_solver.is_infeasible():
            raise ValueError(' μ = 0 is infeasible, check matrix consistency.')
        else:
            raise ValueError(' ' + self.unknown_flag_msg(0))

        # bissection
        mu_min = 0
        mu_max = 2.5
        mu_test = mu_max
        self._sol_basis = None
        iter = 0
        if self.verbose:
            print('  Finding optimal μ at tolerance {} ...'.format(self.bissection_tol))
        while (mu_max - mu_min) > self.bissection_tol and (self.max_bissection_iters is None or iter < self.max_bissection_iters):
            iter += 1

            if self.verbose:
                print('    Iteration {} for μ = {} ...'.format(iter, mu_test), end='')

            self.matrix.build_matrices(mu_test)
            self.lp_solver.build_lp()
            self.lp_solver.solve_lp()
            if self.lp_solver.is_feasible():
                mu_min = mu_test
                self.lp_solver.store_results(mu_test)
                if self.verbose:
                    print(' is feasible.')
            elif self.lp_solver.is_infeasible():
                mu_max = mu_test
                if self.verbose:
                    print(' is infeasible.')
            else:
                raise ValueError(' ' + self.unknown_flag_msg(mu_test))

            # next μ to be tested
            mu_test = (mu_min + mu_max) / 2

        if (mu_max - mu_min) > self.bissection_tol:
            raise ValueError('Optimal solution could not be found in {} iterations\n  mu_min: {}\n  mu_max: {}'.format(
                self.max_bissection_iters, mu_min, mu_max))

        # Recompute the matrices so that it necessarily correspond to the solution
        self.matrix.build_matrices(self.mu_opt)

    def solve_grid(self):
        """Compute configuration corresponding to maximal growth rate."""
        self._sol_basis = None
        self.X = self.lambda_ = self.mu_opt = None

        # check that μ=0 is solution
        self.matrix.build_matrices(0)
        self.lp_solver.build_lp()
        self.lp_solver.solve_lp()
        if self.lp_solver.is_feasible():
            self.lp_solver.store_results(0.0)
        elif self.lp_solver.is_infeasible():
            raise ValueError('μ = 0 is infeasible, check matrix consistency.')
        else:
            raise ValueErro(self.unknown_flag_msg(0))

        # grid
        vec_mu = numpy.arange(0, 0.8, 0.001)
        for mu_test in vec_mu:
            self.matrix.build_matrices(mu_test)
            self.lp_solver.build_lp()
            self.lp_solver.solve_lp()
            if self.lp_solver.is_feasible():
                if self.verbose:
                    print('μ = ' + str(mu_test) + ' is feasible.')
                self.lp_solver.store_results(mu_test)
            elif self.lp_solver.is_infeasible():
                if self.verbose:
                    print('μ = ' + str(mu_test) + ' is infeasible.')
            else:
                raise ValueError(self.unknown_flag_msg(mu_test))

    def unknown_flag_msg(self, mu):
        status = self.lp_solver.get_status()
        return ('At μ = {}: Unknown exit flag {} corresponding to status {}. '
                'Interrupting computation.'.format(
                    mu, status['code'], status['message']
                ))


class LpSolver(abc.ABC):
    """ Base class for LP solver

    Attributes:
        name (:obj:`str`): name
        rba_solver (:obj:`Solver`): RBA solver
    """

    def __init__(self, rba_solver):
        """
        Args:
            rba_solver (:obj:`Solver`): RBA solver
        """
        self.rba_solver = rba_solver

    @property
    @abc.abstractmethod
    def name(self):
        pass

    @abc.abstractmethod
    def build_lp(self):
        """ Build an LP problem based on current matrices. """
        pass

    @abc.abstractmethod
    def solve_lp(self):
        pass

    @abc.abstractmethod
    def is_feasible(self):
        pass

    @abc.abstractmethod
    def is_infeasible(self):
        pass

    @abc.abstractmethod
    def store_results(self, mu):
        pass

    @abc.abstractmethod
    def get_status():
        pass


class CplexLpSolver(LpSolver):
    """ CPLEX LP solver """

    @property
    def name(self):
        return 'cplex'

    def build_lp(self):
        """
        Build CPLEX problem based on current matrices.

        Returns
        -------
        cplex object
            LP problem with default parameters corresponding to current
            matrices.

        """
        import cplex
        rba_solver = self.rba_solver

        # preprocess matrices
        # rescale concentration columns?
        lhs = rba_solver.matrix.A
        # scaling_factor = 1000
        # scaling = numpy.ones(lhs.shape[1])
        # scaling[numpy.concatenate([rba_solver.enzyme_cols, rba_solver.process_cols,
        # rba_solver.species_cols])] \
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
        lp_problem.variables.add(names=rba_solver.matrix.col_names)
        lp_problem.variables.set_lower_bounds(zip(rba_solver.matrix.col_names, rba_solver.matrix.LB))
        lp_problem.variables.set_upper_bounds(zip(rba_solver.matrix.col_names, rba_solver.matrix.UB))
        lp_problem.objective.set_sense(lp_problem.objective.sense.minimize)
        lp_problem.objective.set_linear(zip(rba_solver.matrix.col_names, rba_solver.matrix.f))

        lp_problem.linear_constraints.add(names=rba_solver.matrix.row_names)
        lp_problem.linear_constraints.set_linear_components(zip(rba_solver.matrix.row_names, rows))
        lp_problem.linear_constraints.set_rhs(zip(rba_solver.matrix.row_names, rba_solver.matrix.b))
        lp_problem.linear_constraints.set_senses(zip(rba_solver.matrix.row_names, rba_solver.matrix.row_signs))

        # set starting point (not exactly sure how this works)
        if rba_solver._sol_basis is not None:
            lp_problem.start.set_start(
                rba_solver._sol_basis[0], rba_solver._sol_basis[1], rba_solver.X, [], [],
                rba_solver.lambda_)

        self._lp = lp_problem

    def solve_lp(self):
        self._lp.solve()

    def is_feasible(self):
        return self._lp.solution.get_status() == self._lp.solution.status.optimal

    def is_infeasible(self):
        flag = self._lp.solution.get_status()
        return (flag == self._lp.solution.status.infeasible
                or flag == self._lp.solution.status.optimal_infeasible)

    def get_status(self):
        return {
            'code': self._lp.solution.get_status(),
            'message': self._lp.solution.get_status_string(),
        }

    def store_results(self, mu):
        rba_solver = self.rba_solver
        rba_solver.mu_opt = mu
        rba_solver.X = numpy.array(self._lp.solution.get_values())
        rba_solver.lambda_ = self._lp.solution.get_dual_values()
        rba_solver._sol_basis = self._lp.solution.basis.get_basis()


class OptlangLpSolver(LpSolver):
    """ OptLang LP solver """

    def __init__(self, rba_solver, lp_solver):
        """
        Args:
            rba_solver (:obj:`Solver`): RBA solver
        """
        super(OptlangLpSolver, self).__init__(rba_solver)
        self.lp_solver = lp_solver

    @property
    def name(self):
        return self.lp_solver

    def build_lp(self):
        """ Build an LP problem based on current matrices. """
        import optlang
        solver_interface = getattr(optlang, self.lp_solver + '_interface')

        rba_solver = self.rba_solver

        # define problem
        model = solver_interface.Model()

        # set parameters
        model.configuration.tolerances.feasibility = 1e-9
        model.configuration.tolerances.integrality = 1e-9
        if self.lp_solver == 'cplex':
            model.configuration.tolerances.optimality = 1e-9

        # define columns and add rows
        variables = []
        obj = 0
        for name, lb, ub, obj_coeff in zip(rba_solver.matrix.col_names, rba_solver.matrix.LB,
                                           rba_solver.matrix.UB, rba_solver.matrix.f):
            variable = solver_interface.Variable(name, lb=lb, ub=ub)
            variables.append(variable)
            obj += obj_coeff * variable

        model.objective = solver_interface.Objective(obj, direction='min')

        constraints = rba_solver.matrix.A.tolil()
        for name, i_variables, coeffs, rhs, sense in zip(rba_solver.matrix.row_names, constraints.rows, constraints.data,
                                                         rba_solver.matrix.b, rba_solver.matrix.row_signs):
            constraint = 0

            for i_variable, coeff in zip(i_variables, coeffs):
                constraint += coeff * variables[i_variable]

            if sense == 'E':  # equality
                lb = rhs
                ub = rhs
            elif sense == 'L':  # less than
                lb = None
                ub = rhs
            elif sense == 'G':  # greater than
                lb = rhs
                ub = None
            else:
                raise NotImplementedError('Constraint sense `{}` is not supported.'.format(sense))

            model.add(solver_interface.Constraint(constraint, lb=lb, ub=ub))

        self._model = model

    def solve_lp(self):
        self._model.optimize()

    def is_feasible(self):
        import optlang.interface
        return self._model.status in [optlang.interface.OPTIMAL, optlang.interface.FEASIBLE]

    def is_infeasible(self):
        import optlang.interface
        return self._model.status in [optlang.interface.INFEASIBLE]

    def get_status(self):
        import optlang.interface
        return {
            'code': self._model.status,
            'message': optlang.interface.statuses[self._model.status],
        }

    def store_results(self, mu):
        rba_solver = self.rba_solver
        rba_solver.mu_opt = mu
        rba_solver.X = numpy.array(list(self._model.primal_values.values()))
        rba_solver.lambda_ = list(self._model.shadow_prices.values())
        rba_solver._sol_basis = None
