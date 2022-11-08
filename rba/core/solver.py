"""Module defining Solver class."""

# python 2/3 compatibility
from __future__ import division, print_function, absolute_import

# global imports
import abc
import functools
import numpy
import os
from scipy.sparse import coo_matrix

# from scipy.sparse import diags
import sys
try:
    import conv_opt
except ModuleNotFoundError:
    conv_opt = None
try:
    import optlang
    import optlang.interface
except ModuleNotFoundError:
    optlang = None
try:
    import cplex
except ModuleNotFoundError:
    cplex = None
try:
    import glpk
except ModuleNotFoundError:
    glpk = None
try:
    import gurobipy
except ModuleNotFoundError:
    gurobipy = None
try:
    import swiglpk
except ModuleNotFoundError:
    swiglpk = None


def is_conv_opt_available():
    """ Determine whether ConvOpt is available for solving LP problems

    Returns:
        :obj:`bool`: whether ConvOpt is available for solving LP problems
    """
    return conv_opt is not None


def is_optlang_available():
    """ Determine whether OptLang is available for solving LP problems

    Returns:
        :obj:`bool`: whether OptLang is available for solving LP problems
    """
    return optlang is not None


def is_cplex_available():
    """ Determine whether CPLEX is available for solving LP problems

    Returns:
        :obj:`bool`: whether CPLEX is available for solving LP problems
    """
    return cplex is not None


def is_glpk_available():
    """ Determine whether GLPK is available for solving LP problems

    Returns:
        :obj:`bool`: whether GLPK is available for solving LP problems
    """
    return glpk is not None


def is_swiglpk_available():
    """ Determine whether GLPK is available for solving LP problems

    Returns:
        :obj:`bool`: whether GLPK is available for solving LP problems
    """
    return swiglpk is not None


def get_gurobi_env_vars():
    """ Get environment variables for Gurobi environment variables

    Returns:
        :obj:`dict`: dictionary of Gurobi environment variables
    """
    params = {}
    for key, val in os.environ.items():
        if key.startswith('GRB_') and len(key) > 4:
            key = key[4:]
            if key == 'LICENSEID':
                val = int(val)
            params[key] = val

    return params


def is_gurobi_available():
    """ Determine whether Gurobi is available for solving LP problems

    Returns:
        :obj:`bool`: whether Gurobi is available for solving LP problems
    """
    return (
        gurobipy
        and (
            get_gurobi_env_vars()
            or os.path.isfile(os.path.expanduser(os.path.join('~', 'gurobi.lic')))
            or os.path.isfile(os.path.join('/opt', 'gurobi', 'gurobi.lic'))
        )
    )


class Solver(object):
    """RBA solver."""

    def __init__(self, matrix, lp_solver=None, mu_min=0, mu_max=2.5, bissection_tol=1e-6, max_bissection_iters=None, verbose=False):
        """
        Args:
            matrix
            lp_solver (:obj:`str`, optional): preferred solver
            mu_min (:obj:`float`, optional): minimum μ to check
            mu_max (:obj:`float`, optional): maximum μ to check
            bissection_tol (:obj:`float`, optional): bissection tolerance
            max_bissection_iters (:obj:`int`, optional): maximum number of bissection iterations
            verbose (:obj:`bool`, optional): whether to display diagnostic information
        """
        self.matrix = matrix

        if (
            (is_cplex_available())
                and (lp_solver in ['cplex'] or lp_solver is None)
        ):
            self.lp_solver = CplexLpSolver(self)
        elif (
            (is_cplex_available() and is_conv_opt_available())
            and (lp_solver in ['cplex', 'cplex_conv_opt'] or lp_solver is None)
        ):
            self.lp_solver = ConvOptLpSolver(self, 'cplex')
        elif (
            (is_cplex_available() and is_optlang_available())
            and (lp_solver in ['cplex', 'cplex_optlang'] or lp_solver is None)
        ):
            self.lp_solver = OptlangLpSolver(self, 'cplex')

        elif (
            (is_gurobi_available())
            and (lp_solver in ['gurobi'] or lp_solver is None)
        ):
            self.lp_solver = GurobiLpSolver(self)
        elif (
            (is_gurobi_available() and is_conv_opt_available())
            and (lp_solver in ['gurobi', 'gurobi_conv_opt'] or lp_solver is None)
        ):
            self.lp_solver = ConvOptLpSolver(self, 'gurobi')
        elif (
            (is_gurobi_available() and is_optlang_available())
            and (lp_solver in ['gurobi', 'gurobi_optlang'] or lp_solver is None)
        ):
            self.lp_solver = OptlangLpSolver(self, 'gurobi')
        elif (
            (is_glpk_available())
            and (lp_solver in ['glpk'] or lp_solver is None)
        ):
            self.lp_solver = GlpkLpSolver(self)
        elif (
            (is_swiglpk_available())
            and (lp_solver in ['swiglpk'] or lp_solver is None)
        ):
            self.lp_solver = SwiglpkLpSolver(self)
        elif (
            (is_swiglpk_available() and is_optlang_available())
            and (lp_solver in ['glpk', 'glpk_optlang'] or lp_solver is None)
        ):
            self.lp_solver = OptlangLpSolver(self, 'glpk')
        elif (
            (is_swiglpk_available() and is_conv_opt_available())
            and (lp_solver in ['glpk', 'glpk_conv_opt'] or lp_solver is None)
        ):
            self.lp_solver = ConvOptLpSolver(self, 'glpk')

        elif lp_solver is None or lp_solver == 'scipy':
            self.lp_solver = OptlangLpSolver(self, 'scipy')

        else:
            raise NotImplementedError('LP solver `{}` is not implemented'.format(lp_solver))

        self.mu_min = mu_min
        self.mu_max = mu_max
        self.bissection_tol = bissection_tol
        self.max_bissection_iters = max_bissection_iters

        self.verbose = verbose

    def solve(self):
        """ Compute configuration corresponding to maximal growth rate. """
        self._sol_basis = None
        self.X = self.lambda_ = self.mu_opt = None
        mu_min = self.mu_min
        mu_max = self.mu_max

        # check that μ=μ_min is solution
        if self.verbose:
            print('  Checking μ = μ_min = {} is feasible ...'.format(mu_min), end='')
            sys.stdout.flush()
        self.matrix.build_matrices(mu_min)
        self.lp_solver.build_lp()
        self.lp_solver.solve_lp()
        if self.lp_solver.is_feasible():
            self.lp_solver.store_results(mu_min)
            if self.verbose:
                print(' μ = μ_min = {} is feasible.'.format(mu_min))
        elif self.lp_solver.is_infeasible():
            raise ValueError(' μ = μ_min = {} is infeasible, check matrix consistency.'.format(mu_min))
        else:
            raise ValueError(' ' + self.unknown_flag_msg(mu_min))

        # bissection
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
                    print(' is feasible. At status {} ({})'.format(self.lp_solver.get_status()['code'],self.lp_solver.get_status()['message']))
            elif self.lp_solver.is_infeasible():
                mu_max = mu_test
                if self.verbose:
                    print(' is infeasible. At status {} ({})'.format(self.lp_solver.get_status()['code'],self.lp_solver.get_status()['message']))
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
        """ Compute configuration corresponding to maximal growth rate. """
        self._sol_basis = None
        self.X = self.lambda_ = self.mu_opt = None
        mu_min = self.mu_min
        mu_max = self.mu_max

        # check that μ=0 is solution
        self.matrix.build_matrices(mu_min)
        self.lp_solver.build_lp()
        self.lp_solver.solve_lp()
        if self.lp_solver.is_feasible():
            self.lp_solver.store_results(mu_min)
        elif self.lp_solver.is_infeasible():
            raise ValueError('μ = μ_min = {} is infeasible, check matrix consistency.'.format(mu_min))
        else:
            raise ValueError(self.unknown_flag_msg(mu_min))

        # grid
        # TODO: parallelize with multiprocessing module
        vec_mu = numpy.arange(mu_min, mu_max, 0.001)
        for mu_test in vec_mu:
            self.matrix.build_matrices(mu_test)
            self.lp_solver.build_lp()
            self.lp_solver.solve_lp()
            if self.lp_solver.is_feasible():
                if self.verbose:
                    print('μ = {} is feasible. At status {} ({})'.format(str(mu_test),self.lp_solver.get_status()['code'],self.lp_solver.get_status()['message']))
                self.lp_solver.store_results(mu_test)
            elif self.lp_solver.is_infeasible():
                if self.verbose:
                    print('μ = {} is infeasible. At status {} ({})'.format(str(mu_test),self.lp_solver.get_status()['code'],self.lp_solver.get_status()['message']))
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

    def init_solver(self, lp_solver):
        if lp_solver == 'gurobi' and not getattr(gurobipy, '__initialized__', None):
            params = get_gurobi_env_vars()
            if params:
                env = gurobipy.Env(params=params)
                gurobipy.Model = functools.partial(gurobipy.Model, env=env)  # to make OptLang use the Gurobi environment
                gurobipy.__initialized__ = True

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
        #import cplex
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


class SwiglpkLpSolver(LpSolver):
    """ Swiglpk LP solver

    Attributes:
        rba_solver (:obj:`Solver`): RBA solver
        verbose (:obj:`bool`): whether to display diagnostic information
    """

    def __init__(self, rba_solver, verbose=False):
        """
        Args:
            rba_solver (:obj:`Solver`): RBA solver
            verbose (:obj:`bool`, optional): whether to display diagnostic information
        """
        super(SwiglpkLpSolver, self).__init__(rba_solver)
        self.verbose = verbose

    @property
    def name(self):
        return 'swiglpk'

    def build_lp(self):
        """
        Constructs a GLPK-object.
        """
        rba_solver = self.rba_solver

        # define problem
        swiglpk.glp_term_out(swiglpk.GLP_OFF)
        lp = swiglpk.glp_create_prob()
        swiglpk.glp_create_index(lp)
        swiglpk.glp_scale_prob(lp, swiglpk.GLP_SF_AUTO)

        swiglpk.glp_set_obj_dir(lp, swiglpk.GLP_MIN)
        swiglpk.glp_add_rows(lp, len(rba_solver.matrix.row_names)) #  3=nrs number of constraints
        swiglpk.glp_add_cols(lp, len(rba_solver.matrix.col_names)) #  3=nrs number of variables

        row_sign_mapping={'E':swiglpk.GLP_FX,'L':swiglpk.GLP_UP,'G':swiglpk.GLP_LO}
        for row_index in range(len(rba_solver.matrix.row_names)):
            swiglpk.glp_set_row_name(lp, row_index+1, rba_solver.matrix.row_names[row_index]) #sets name of first row to p -> note that indexing starts at 1 not 0
            swiglpk.glp_set_row_bnds(lp, row_index+1, row_sign_mapping[rba_solver.matrix.row_signs[row_index]], rba_solver.matrix.b[row_index], rba_solver.matrix.b[row_index]) #bounds the first row between lb 0 and ub 100

        bound_type_map={False:swiglpk.GLP_DB,True:swiglpk.GLP_FX}
        for col_index in range(len(rba_solver.matrix.col_names)):
            swiglpk.glp_set_col_name(lp, col_index+1, rba_solver.matrix.col_names[col_index])
            swiglpk.glp_set_obj_coef(lp, col_index+1, float(rba_solver.matrix.f[col_index])) #sets name of first variables objective coefficient to to 10 -> note that indexing starts at 1 not 0
            lb=rba_solver.matrix.LB[col_index]
            ub=rba_solver.matrix.UB[col_index]
            swiglpk.glp_set_col_bnds(lp, col_index+1, bound_type_map[lb==ub], float(lb), float(ub))

        constraints = rba_solver.matrix.A.tocoo()
        nonzero_row_inds=list(constraints.row)
        nonzero_col_inds=list(constraints.col)
        nonzero_data=list(constraints.data)

        nonzero_row_indices = swiglpk.intArray(len(nonzero_data)+1)
        nonzero_col_indices = swiglpk.intArray(len(nonzero_data)+1)
        nonzero_coefficients = swiglpk.doubleArray(len(nonzero_data)+1)

        for i in range(len(nonzero_data)):
            nonzero_row_indices[i+1]=int(nonzero_row_inds[i]+1)
            nonzero_col_indices[i+1]=int(nonzero_col_inds[i]+1)
            nonzero_coefficients[i+1]=nonzero_data[i]

        swiglpk.glp_load_matrix(lp, len(nonzero_data), nonzero_row_indices, nonzero_col_indices, nonzero_coefficients)

        swiglpk.glp_scale_prob(lp,swiglpk.GLP_SF_EQ) #GLP_SF_EQ
        #GLP_SF_GM perform geometric mean scaling;
        #GLP_SF_EQ perform equilibration scaling;
        #GLP_SF_2N round scale factors to nearest power of two;
        #GLP_SF_SKIP skip scaling, if the problem is well scaled.

        #Defining GLPK solver parameters
        self.glpk_simplex_params=swiglpk.glp_smcp()
        setattr(self.glpk_simplex_params, "tol_bnd", 1e-9)
        setattr(self.glpk_simplex_params, "tol_dj", 1e-9)
        setattr(self.glpk_simplex_params, "tol_piv", 1e-9)
        #setattr(self.glpk_simplex_params, "presolve", swiglpk.GLP_OFF)
        #setattr(self.glpk_simplex_params, "tm_lim", 2000)

        swiglpk.glp_init_smcp(self.glpk_simplex_params)
        self._model = lp

    def solve_lp(self):
        swiglpk.glp_simplex(self._model, self.glpk_simplex_params)
        #swiglpk.glp_simplex(self._model,None)

    def is_feasible(self):
        #return swiglpk.glp_get_status(self._model) in [swiglpk.GLP_FEAS,swiglpk.GLP_OPT,swiglpk.GLP_UNDEF]
        return swiglpk.glp_get_status(self._model) in [swiglpk.GLP_FEAS,swiglpk.GLP_OPT]

    def is_infeasible(self):
        #return swiglpk.glp_get_status(self._model) not in [swiglpk.GLP_FEAS,swiglpk.GLP_OPT,swiglpk.GLP_UNDEF]
        return swiglpk.glp_get_status(self._model) not in [swiglpk.GLP_FEAS,swiglpk.GLP_OPT]

    def get_status(self):
        return {
            'code': swiglpk.glp_get_status(self._model),
            'message': self.get_status_message(swiglpk.glp_get_status(self._model)),
        }

    def get_status_message(self, code):
        codes = {
            swiglpk.GLP_OPT: 'solution is optimal',
            swiglpk.GLP_UNDEF: 'solution is undefined (infeasibility is assumed)',
            swiglpk.GLP_FEAS: 'solution is feasible, but not necessarily optimal',
            swiglpk.GLP_INFEAS: 'solution is infeasible',
            swiglpk.GLP_NOFEAS: 'problem has no feasible solution',
            swiglpk.GLP_UNBND: 'problem has an unbounded solution',
        }
        return codes.get(code, None)

    def store_results(self, mu):
        rba_solver = self.rba_solver
        rba_solver.mu_opt = mu
        rba_solver.X = numpy.array([swiglpk.glp_get_col_prim(self._model, i+1) for i in range(len(rba_solver.matrix.col_names))])
        rba_solver.lambda_ = [swiglpk.glp_get_row_dual(self._model, i+1) for i in range(len(rba_solver.matrix.row_names))]
        rba_solver._sol_basis = None


class GlpkLpSolver(LpSolver):
    """ GLPK LP solver

    Attributes:
        rba_solver (:obj:`Solver`): RBA solver
        verbose (:obj:`bool`): whether to display diagnostic information
    """

    def __init__(self, rba_solver, verbose=False):
        """
        Args:
            rba_solver (:obj:`Solver`): RBA solver
            verbose (:obj:`bool`, optional): whether to display diagnostic information
        """
        super(GlpkLpSolver, self).__init__(rba_solver)
        self.verbose = verbose

    @property
    def name(self):
        return 'glpk'

    def build_lp(self):
        """ Build an LP problem based on current matrices. """
        rba_solver = self.rba_solver

        # define problem
        lp = glpk.LPX()
        lp.rows.add(len(rba_solver.matrix.row_names))
        lp.cols.add(len(rba_solver.matrix.col_names))

        # define variables and objective
        for i_variable, (lb, ub, obj_coeff) in enumerate(zip(rba_solver.matrix.LB, rba_solver.matrix.UB, rba_solver.matrix.f)):
            lp.cols[i_variable].bounds = lb, ub
            lp.obj[i_variable] = obj_coeff

        lp.obj.maximize = False

        # define constraints
        constraints = rba_solver.matrix.A.tolil()
        matrix = []
        for i_constraint, (i_variables, coeffs, rhs, sense) in enumerate(zip(constraints.rows, constraints.data,
                                                                             rba_solver.matrix.b, rba_solver.matrix.row_signs)):

            for i_variable, coeff in zip(i_variables, coeffs):
                matrix.append((i_constraint, i_variable, coeff))

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

            lp.rows[i_constraint].bounds = lb, ub

        lp.matrix = matrix

        # define options
        lp.scale(glpk.LPX.SF_GM)
        if self.verbose:
            glpk.env.term_hook = None
        else:
            glpk.env.term_hook = lambda output: None

        self._model = lp

    def solve_lp(self):
        self._model.simplex(tol_bnd=1e-9, tol_dj=1e-9, tol_piv=1e-10)

    def is_feasible(self):
        return self._model.status == 'opt'

    def is_infeasible(self):
        return self._model.status in ['infeas', 'nofeas']

    def get_status(self):
        return {
            'code': self._model.status,
            'message': self.get_status_message(self._model.status),
        }

    def get_status_message(self, code):
        codes = {
            'opt': 'solution is optimal',
            'undef': 'solution is undefined',
            'feas': 'solution is feasible, but not necessarily optimal',
            'infeas': 'solution is infeasible',
            'nofeas': 'problem has no feasible solution',
            'unbnd': 'problem has an unbounded solution',
        }
        return codes.get(code, None)

    def store_results(self, mu):
        rba_solver = self.rba_solver
        rba_solver.mu_opt = mu
        rba_solver.X = numpy.array([row.primal_s for row in self._model.rows])
        rba_solver.lambda_ = [col.dual_s for col in self._model.cols]
        rba_solver._sol_basis = None


class GurobiLpSolver(LpSolver):
    """ Gurobi LP solver

    Attributes:
        rba_solver (:obj:`Solver`): RBA solver
        verbose (:obj:`bool`): whether to display diagnostic information
    """

    def __init__(self, rba_solver, verbose=False):
        """
        Args:
            rba_solver (:obj:`Solver`): RBA solver
            verbose (:obj:`bool`, optional): whether to display diagnostic information
        """
        super(GurobiLpSolver, self).__init__(rba_solver)
        self.verbose = verbose

        self.init_solver('gurobi')

    @property
    def name(self):
        return 'gurobi'

    def build_lp(self):
        """ Build an LP problem based on current matrices. """
        rba_solver = self.rba_solver

        # define problem
        model = gurobipy.Model()

        # define variables
        variables = model.addVars(len(rba_solver.matrix.col_names),
                                  lb=rba_solver.matrix.LB,
                                  ub=rba_solver.matrix.UB,
                                  obj=rba_solver.matrix.f,
                                  vtype=gurobipy.GRB.CONTINUOUS,
                                  )

        model.setObjective(model.getObjective(), gurobipy.GRB.MINIMIZE)

        # define constraints
        constraints = rba_solver.matrix.A.tolil()

        for i_variables, coeffs, rhs, sense in zip(constraints.rows, constraints.data,
                                                   rba_solver.matrix.b, rba_solver.matrix.row_signs):

            vars = []
            for i_variable, coeff in zip(i_variables, coeffs):
                vars.append(variables[i_variable])

            if sense == 'E':  # equality
                sense = gurobipy.GRB.EQUAL
            elif sense == 'L':  # less than
                sense = gurobipy.GRB.LESS_EQUAL
            elif sense == 'G':  # greater than
                sense = gurobipy.GRB.GREATER_EQUAL
            else:
                raise NotImplementedError('Constraint sense `{}` is not supported.'.format(sense))

            model.addLConstr(gurobipy.LinExpr(coeffs, vars), sense, rhs=rhs)

        # define options
        model.setParam('FeasibilityTol', 1e-9)
        model.setParam('OptimalityTol', 1e-9)
        model.setParam('MarkowitzTol', 0.1)
        model.setParam("ScaleFlag", 3)
        if self.verbose:
            model.setParam('LogToConsole', 1)
        else:
            model.setParam('LogToConsole', 0)

        self._model = model

    def solve_lp(self):
        self._model.optimize()

    def is_feasible(self):
        return self._model.status == gurobipy.GRB.Status.OPTIMAL

    def is_infeasible(self):
        return self._model.status == gurobipy.GRB.Status.INFEASIBLE

    def get_status(self):
        return {
            'code': self._model.status,
            'message': self.get_status_message(self._model.status),
        }

    def get_status_message(self, code):
        codes = {
            gurobipy.GRB.LOADED: (
                'Model is loaded, but no solution information is available.'
            ),
            gurobipy.GRB.OPTIMAL: (
                'Model was solved to optimality (subject to tolerances), and an optimal solution is available.'
            ),
            gurobipy.GRB.INFEASIBLE: (
                'Model was proven to be infeasible.'
            ),
            gurobipy.GRB.INF_OR_UNBD: (
                'Model was proven to be either infeasible or unbounded. To obtain a more definitive conclusion, '
                'set the DualReductions parameter to 0 and reoptimize.'
            ),
            gurobipy.GRB.UNBOUNDED: (
                'Model was proven to be unbounded. Important note: an unbounded status indicates the presence of an unbounded ray '
                'that allows the objective to improve without limit. It says nothing about whether the model has a feasible solution. '
                'If you require information on feasibility, you should set the objective to zero and reoptimize.'
            ),
            gurobipy.GRB.CUTOFF: (
                'Optimal objective for model was proven to be worse than the value specified in the Cutoff parameter. '
                'No solution information is available.'
            ),
            gurobipy.GRB.ITERATION_LIMIT: (
                'Optimization terminated because the total number of simplex iterations performed exceeded the value specified in the '
                'IterationLimit parameter, or because the total number of barrier iterations exceeded the value specified in the '
                'BarIterLimit parameter.'
            ),
            gurobipy.GRB.NODE_LIMIT: (
                'Optimization terminated because the total number of branch-and-cut nodes explored exceeded the value specified '
                'in the NodeLimit parameter.'
            ),
            gurobipy.GRB.TIME_LIMIT: (
                'Optimization terminated because the time expended exceeded the value specified in the TimeLimit parameter.'
            ),
            gurobipy.GRB.SOLUTION_LIMIT: (
                'Optimization terminated because the number of solutions found reached the value specified in the SolutionLimit parameter.'
            ),
            gurobipy.GRB.INTERRUPTED: (
                'Optimization was terminated by the user.'
            ),
            gurobipy.GRB.NUMERIC: (
                'Optimization was terminated due to unrecoverable numerical difficulties.'
            ),
            gurobipy.GRB.SUBOPTIMAL: (
                'Unable to satisfy optimality tolerances; a sub-optimal solution is available.'
            ),
            gurobipy.GRB.INPROGRESS: (
                'An asynchronous optimization call was made, but the associated optimization run is not yet complete.'
            ),
            gurobipy.GRB.USER_OBJ_LIMIT: (
                'User specified an objective limit (a bound on either the best objective or the best bound), and that '
                'limit has been reached.'
            ),
        }
        return codes.get(code, None)

    def store_results(self, mu):
        rba_solver = self.rba_solver
        rba_solver.mu_opt = mu
        rba_solver.X = numpy.array([var.getAttr('X') for var in self._model.getVars()])
        rba_solver.lambda_ = [constraint.getAttr('Pi') for constraint in self._model.getConstrs()]
        rba_solver._sol_basis = None


class ConvOptLpSolver(LpSolver):
    """ ConvOpt LP solver

    Attributes:
        rba_solver (:obj:`Solver`): RBA solver
        lp_solver (:obj:`str`): LP solver
        verbose (:obj:`bool`): whether to display diagnostic information
    """

    def __init__(self, rba_solver, lp_solver, verbose=False):
        """
        Args:
            rba_solver (:obj:`Solver`): RBA solver
            lp_solver (:obj:`str`): LP solver
            verbose (:obj:`bool`, optional): whether to display diagnostic information
        """
        super(ConvOptLpSolver, self).__init__(rba_solver)
        self.lp_solver = lp_solver
        self.verbose = verbose

        self.init_solver(lp_solver)

    @property
    def name(self):
        return self.lp_solver

    def build_lp(self):
        """ Build an LP problem based on current matrices. """
        rba_solver = self.rba_solver

        # define problem
        model = conv_opt.Model()

        # define columns and add rows
        for name, lb, ub, obj_coeff in zip(rba_solver.matrix.col_names, rba_solver.matrix.LB,
                                           rba_solver.matrix.UB, rba_solver.matrix.f):
            var = conv_opt.Variable(name=name, type=conv_opt.VariableType.continuous, lower_bound=lb, upper_bound=ub)
            model.variables.append(var)

            model.objective_terms.append(conv_opt.LinearTerm(var, obj_coeff))

        model.objective_direction = conv_opt.ObjectiveDirection.minimize

        constraints = rba_solver.matrix.A.tolil()
        for name, i_variables, coeffs, rhs, sense in zip(rba_solver.matrix.row_names, constraints.rows, constraints.data,
                                                         rba_solver.matrix.b, rba_solver.matrix.row_signs):

            terms = []
            for i_variable, coeff in zip(i_variables, coeffs):
                var = model.variables[i_variable]
                terms.append(conv_opt.LinearTerm(var, coeff))

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

            constraint = conv_opt.Constraint(terms, name=name, lower_bound=lb, upper_bound=ub)
            model.constraints.append(constraint)

        self._model = model

    def solve_lp(self):
        # set parameters
        self._options = conv_opt.SolveOptions(
            solver=getattr(conv_opt.Solver, self.lp_solver),
            presolve=conv_opt.Presolve.off,
            verbosity=conv_opt.Verbosity.status if self.verbose else conv_opt.Verbosity.off,
        )

        self._result = self._model.solve(self._options)

    def is_feasible(self):
        return self._result.status_code in [conv_opt.StatusCode.optimal]

    def is_infeasible(self):
        return self._result.status_code in [conv_opt.StatusCode.infeasible]

    def get_status(self):
        return {
            'code': self._result.status_code,
            'message': self._result.status_message,
        }

    def store_results(self, mu):
        rba_solver = self.rba_solver
        rba_solver.mu_opt = mu
        rba_solver.X = self._result.primals
        rba_solver.lambda_ = self._result.duals.tolist()
        rba_solver._sol_basis = None


class OptlangLpSolver(LpSolver):
    """ OptLang LP solver

    Attributes:
        rba_solver (:obj:`Solver`): RBA solver
        lp_solver (:obj:`str`): LP solver
        verbose (:obj:`bool`): whether to display diagnostic information
    """

    def __init__(self, rba_solver, lp_solver, verbose=False):
        """
        Args:
            rba_solver (:obj:`Solver`): RBA solver
            lp_solver (:obj:`str`): LP solver
            verbose (:obj:`bool`, optional): whether to display diagnostic information
        """
        super(OptlangLpSolver, self).__init__(rba_solver)
        self.lp_solver = lp_solver
        self.verbose = verbose

        self.init_solver(lp_solver)

    @property
    def name(self):
        return self.lp_solver

    def build_lp(self):
        """ Build an LP problem based on current matrices. """
        solver_interface = getattr(optlang, self.lp_solver + '_interface')

        rba_solver = self.rba_solver

        # define problem
        model = solver_interface.Model()

        # set parameters
        if self.lp_solver in ['cplex']:
            model.configuration.tolerances.feasibility = 1e-9
            model.configuration.tolerances.integrality = 1e-9
            model.configuration.tolerances.optimality = 1e-9
        elif self.lp_solver in ['glpk']:
            model.configuration.tolerances.feasibility = 1e-9
            model.configuration.tolerances.integrality = 1e-9
        elif self.lp_solver in ['gurobi']:
            model.configuration.tolerances.feasibility = 1e-9
            model.configuration.tolerances.integrality = 1e-9
            model.configuration.tolerances.optimality = 1e-9

        model.configuration.verbosity = 1 if self.verbose else 0

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
        return self._model.status in [optlang.interface.OPTIMAL, optlang.interface.FEASIBLE]

    def is_infeasible(self):
        return self._model.status in [optlang.interface.INFEASIBLE]

    def get_status(self):
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
