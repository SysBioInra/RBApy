import cplex
import numpy
from scipy.sparse import coo_matrix, diags, hstack, vstack, eye

class RbaSolver(object):
    def __init__(self, blocks):
        self._blocks = blocks
        # convenience variables
        reactions = blocks.reactions
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
        self.col_names = reactions + enzymes \
                         + [p + '_machinery' for p in processes] \
                         + [m + '_flux' for m in undetermined_fluxes]
        self.reaction_cols = numpy.arange(nb_reactions)
        self.enzyme_cols = self.reaction_cols[-1]+1+numpy.arange(nb_enzymes)
        self.process_cols = self.enzyme_cols[-1]+1+numpy.arange(nb_processes)
        self.species_cols = self.process_cols[-1]+1+numpy.arange(nb_undetermined)
        # row information
        self.row_names = blocks.metabolites \
                         + [p + '_capacity' for p in processes] \
                         + [e + '_forward_capacity' for e in enzymes] \
                         + [e + '_backward_capacity' for e in enzymes]\
                         + [c + '_density' for c in compartments]
        self.row_signs = ['E'] * len(blocks.metabolites) \
                         + blocks.processes.capacity_signs \
                         + ['L'] * 2 * nb_enzymes \
                         + blocks.density.signs
        # upper bound, lower bound, objective function
        nb_cols = len(self.col_names)
        self.UB = 1e5 * numpy.ones(nb_cols)
        self.LB = numpy.zeros(nb_cols)
        self.f = numpy.zeros(nb_cols)
        self.LB[self.reaction_cols] = -1e3*numpy.array(blocks.reversibility)
        self.f[self.enzyme_cols] = 1
        # constant building blocks
        self._empty_ExPU = coo_matrix((nb_enzymes, nb_processes+nb_undetermined))
        self._empty_PxR = coo_matrix((nb_processes, nb_reactions))
        self._empty_CxR = coo_matrix((nb_compartments, nb_reactions))
        self._empty_2E = numpy.zeros(2*nb_enzymes)
        # indicator matrices
        R_ind = [reactions.index(r) for r in blocks.enzymes.reaction_catalyzed]
        self._R_to_E = coo_matrix(([1]*nb_enzymes, (range(nb_enzymes), R_ind)),
                                  shape = (nb_enzymes, nb_reactions))
        self._target_reactions \
            = [reactions.index(r) for r in blocks.processes.target_reaction]
        self._lb_reactions \
            = [reactions.index(r) for r in blocks.processes.lb_reaction]
        self._ub_reactions \
            = [reactions.index(r) for r in blocks.processes.ub_reaction]

    def build_matrices(self, mu):        
        ## build A
        # mu-dependent data
        u_composition, u_proc_cost, u_weight \
            = self._blocks.processes.undetermined_values.matrices(mu)
        capacity = self._blocks.processes.capacity.compute(mu)
        (forward, backward) = self._blocks.enzymes.efficiency.compute(mu)
        # stoichiometry constraints
        metab_rows = hstack([self._blocks.S,
                             mu*self._blocks.enzymes.machinery.composition,
                             mu*self._blocks.processes.machinery.composition,
                             u_composition])
        # capacity constraints
        process_rows = hstack([self._empty_PxR,
                               mu*self._blocks.enzymes.machinery.processing_cost,
                               mu*self._blocks.processes.machinery.processing_cost
                               -diags(capacity),
                               u_proc_cost])
        forward_rows = hstack([self._R_to_E, -diags(forward), self._empty_ExPU])
        backward_rows = hstack([-self._R_to_E,
                                -diags(backward), self._empty_ExPU])
        # density constraints
        c_indices = self._blocks.density.compartment_indices
        density_rows = hstack([self._empty_CxR,
                               self._blocks.enzymes.machinery.weight[c_indices],
                               self._blocks.processes.machinery.weight[c_indices],
                               u_weight[c_indices]])
        self.A = vstack([metab_rows, process_rows,
                         forward_rows, backward_rows, density_rows])

        ## build b
        # gather mu-dependent data
        (fluxes, processing, weight) \
            = self._blocks.processes.target_values.compute(mu)
        density_rows = self._blocks.density.values.compute(mu) \
                       - weight[c_indices].T
        # build vector
        self.b = numpy.concatenate([-fluxes, -processing,
                                    self._empty_2E, density_rows])

        ## update lower bounds and upper bounds
        # undetermined metabolites
        self.LB[self.species_cols] \
            = self._blocks.processes.undetermined_values.lb(mu)
        self.UB[self.species_cols] \
            = self._blocks.processes.undetermined_values.ub(mu)
        # target reactions
        self.LB[self.reaction_cols[self._lb_reactions]] \
            = self._blocks.processes.lb.compute(mu)
        self.UB[self.reaction_cols[self._ub_reactions]] \
            = self._blocks.processes.ub.compute(mu)
        r_fluxes = self._blocks.processes.reaction_value.compute(mu)
        self.LB[self.reaction_cols[self._target_reactions]] = r_fluxes
        self.UB[self.reaction_cols[self._target_reactions]] = r_fluxes

    def solve(self):
        """
        Find maximal growth rate.
        """ 
        self._sol_basis = None

        # check that mu=0 is solution
        self.build_matrices(0)
        lp = self.setup_lp()
        lp.solve()
        exit_flag = lp.solution.get_status()
        if (exit_flag == lp.solution.status.optimal):
            self.X = numpy.array(lp.solution.get_values())
            self.lambda_ = lp.solution.get_dual_values()
            self.mu_opt = 0
            self._sol_basis = lp.solution.basis.get_basis()
        elif exit_flag == lp.solution.status.infeasible \
             or (exit_flag == lp.solution.status.optimal_infeasible):
            print('Mu = 0 is infeasible, check matrix consistency')
            return
        else:
            print('At mu = 0: Unknown exit flag '
                  + str(exit_flag) + ' corresponding to status '
                  + lp.solution.get_status_string() + '. '
                  'Interrupting execution')
            return

        # bissection
        mu_min = 0
        mu_max = 2.5
        mu_test = mu_max
        while (mu_max - mu_min) > 1e-4:
            self.build_matrices(mu_test)
            lp = self.setup_lp()
            lp.solve()
            exit_flag = lp.solution.get_status()
            print(mu_test, lp.solution.get_status_string())
            if (exit_flag == lp.solution.status.optimal):
                mu_min = mu_test
                self.X = numpy.array(lp.solution.get_values())
                self.lambda_ = lp.solution.get_dual_values()
                self.mu_opt = mu_min
                self._sol_basis = lp.solution.basis.get_basis()
            elif exit_flag == lp.solution.status.infeasible \
                 or (exit_flag == lp.solution.status.optimal_infeasible):
                mu_max = mu_test
            else:
                print('At mu = ' + str(mu_test) + ': Unknown exit flag '
                      + str(exit_flag) + ' corresponding to status '
                      + lp.solution.get_status_string() + '. '
                      'Interrupting execution')
                return
            # next mu to be tested
            mu_test = (mu_min+mu_max)/2
        print(mu_min)

    def setup_lp(self):
        """
        Setup CPLEX problem based on current matrices.
        """        
        ## preprocess matrices
        # rescale concentration columns?
        lhs = self.A
        #scaling_factor = 1000
        #scaling = numpy.ones(lhs.shape[1])
        #scaling[numpy.concatenate([self.enzyme_cols, self.process_cols,
        #self.species_cols])] \
        #= 1.0/scaling_factor
        #lhs *= diags(scaling)
        
        # transform inequality and equality constraints in CPLEX row format
        lhs = lhs.tolil()
        rows = []
        for nz_ind, data in zip(lhs.rows, lhs.data):
            rows.append(cplex.SparsePair(nz_ind, data))

        ## define problem
        lp_problem = cplex.Cplex()
        # set parameters
        lp_problem.objective.set_sense(lp_problem.objective.sense.minimize)
        lp_problem.parameters.feasopt.tolerance.set(1e-9);
        lp_problem.parameters.simplex.tolerances.feasibility.set(1e-9);
        lp_problem.parameters.simplex.tolerances.optimality.set(1e-9);
        lp_problem.parameters.simplex.tolerances.markowitz.set(0.1);
        lp_problem.parameters.barrier.convergetol.set(1e-9);
        # agressive scaling
        lp_problem.parameters.read.scale.set(1)
        # Threads: the default (0) means that Cplex decides automatically
        # how many threads to use
        # lp_problem.parameters.threads.set(0)
        lp_problem.set_results_stream(None)
        # define columns and add rows
        lp_problem.variables.add(obj = self.f, ub = self.UB, lb = self.LB,
                                 names = self.col_names)
        lp_problem.linear_constraints.add(lin_expr = rows,
                                          rhs = self.b,
                                          senses = self.row_signs,
                                          names = self.row_names)
        # set starting point (not sure how this works)
        if self._sol_basis is not None:
            lp_problem.start.set_basis(self._sol_basis[0], self._sol_basis[1])
            #lp_problem.start.set_start(self._sol_basis[0], self._sol_basis[1],
            #self.X, [], [], self.lambda_)
                                       
        return lp_problem

