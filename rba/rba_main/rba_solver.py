import cplex
import numpy
from scipy.sparse import coo_matrix, diags, hstack, vstack, eye

class RbaSolver(object):
    def __init__(self, model):
        self._model = model
        # convenience variables
        reactions = model.reactions
        enzymes = model.enzymes.ids
        processes = model.processes.ids
        compartments = model.density.compartments
        nb_reactions = len(reactions)
        nb_enzymes = len(enzymes)
        nb_processes = len(processes)
        nb_compartments = len(compartments)
        # column information
        self.col_names = reactions + enzymes \
                         + [p + '_machinery' for p in processes]
        self.reaction_cols = numpy.arange(nb_reactions)
        self.enzyme_cols = nb_reactions + numpy.arange(nb_enzymes)
        self.process_cols = nb_reactions+nb_enzymes+numpy.arange(nb_processes)
        # equality row information
        self.eq_row_names = model.metabolites \
                            + [p + '_capacity' for p in processes] \
                            + [n + '_flux' for n in model.processes.target_reaction]
        # inequality row information
        self.ineq_row_names = [e + '_forward_capacity' for e in enzymes] \
                              + [e + '_backward_capacity' for e in enzymes]\
                              + [c + '_density' for c in compartments]
        # upper bound, lower bound, objective function
        nb_cols = len(self.col_names)
        self.UB = 1e5 * numpy.ones(nb_cols)
        self.LB = numpy.zeros(nb_cols)
        self.f = numpy.zeros(nb_cols)
        self.LB[self.reaction_cols] = -1e3*numpy.array(model.reversibility)
        self.f[self.enzyme_cols] = 1
        # constant building blocks
        self._empty_ExP = coo_matrix((nb_enzymes, nb_processes))
        self._empty_PxR = coo_matrix((nb_processes, nb_reactions))
        self._empty_CxR = coo_matrix((nb_compartments, nb_reactions))
        self._empty_2E = numpy.zeros(2*nb_enzymes)
        # indicator matrices
        R_ind = [reactions.index(r) for r in model.enzymes.reaction_catalyzed]
        self._R_to_E = coo_matrix(([1]*nb_enzymes, (range(nb_enzymes), R_ind)),
                                  shape = (nb_enzymes, nb_reactions))
        R_ind = [reactions.index(r) for r in model.processes.target_reaction]
        nb_TR = len(R_ind)
        self._target_reactions = coo_matrix(([1]*nb_TR, (range(nb_TR), R_ind)),
                                            shape = (nb_TR, nb_cols))
        self._lb_reactions \
            = [reactions.index(r) for r in model.processes.lb_reaction]
        self._ub_reactions \
            = [reactions.index(r) for r in model.processes.ub_reaction]

    def build_matrices(self, mu):
        ## build A
        (forward, backward) = self._model.enzymes.efficiency.compute(mu)
        forward_rows = hstack([self._R_to_E, -diags(forward), self._empty_ExP])
        backward_rows = hstack([-self._R_to_E,
                                -diags(backward), self._empty_ExP])
        c_indices = self._model.density.compartment_indices
        density_rows = hstack([self._empty_CxR,
                               self._model.enzymes.machinery.weight[c_indices],
                               self._model.processes.machinery.weight[c_indices]])
        self.A = vstack([forward_rows, backward_rows, density_rows])
        
        ## build Aeq
        metab_rows = hstack([self._model.S,
                             mu*self._model.enzymes.machinery.composition,
                             mu*self._model.processes.machinery.composition])
        capacity = self._model.processes.capacity.compute(mu)
        process_rows = hstack([self._empty_PxR,
                               mu*self._model.enzymes.machinery.processing_cost,
                               mu*self._model.processes.machinery.processing_cost
                               -diags(capacity)])
        self.Aeq = vstack([metab_rows, process_rows, self._target_reactions])

        ## build b and beq
        (fluxes, processing, weight) \
            = self._model.processes.target_values.compute(mu)
        density_rows = self._model.density.maximum.compute(mu) \
                       - weight[c_indices].T
        r_fluxes = self._model.processes.reaction_value.compute(mu)
        self.b = numpy.concatenate([self._empty_2E, density_rows])
        self.beq = numpy.concatenate([-fluxes, -processing, r_fluxes])

        ## update lower bounds and upper bounds
        self.LB[self.reaction_cols[self._lb_reactions]] \
            = self._model.processes.lb.compute(mu)
        self.UB[self.reaction_cols[self._ub_reactions]] \
            = self._model.processes.ub.compute(mu)

    def solve(self, scaling_factor = 1000):
        """
        Solve the RBA model.
        :param scaling_factor: scaling factor for concentrations vs fluxes
        (used for matrix reconditionning).
        """        
        # cplex exit flags
        OPTIMAL = 1
        INFEASIBLE = 3
        OPTIMAL_UNSCALED_INFEASIBILITIES = 5
        
        # check that mu=0 is solution
        self.build_matrices(0)
        lp_problem = self._setup_lp(scaling_factor)
        lp_problem.solve()
        if lp_problem.solution.get_status() == INFEASIBLE:
            print('Mu = 0 is infeasible, check matrix consistency')
            return self

        # bissection
        mu_min = 0
        mu_max = 2.5
        mu_test = mu_max
        while (mu_max - mu_min) > 1e-4:
            self.build_matrices(mu_test)
            lp_problem = self._setup_lp(scaling_factor)
            lp_problem.solve()
            exit_flag = lp_problem.solution.get_status()
            print(mu_test, lp_problem.solution.get_status_string())
            if (exit_flag == OPTIMAL):
                mu_min = mu_test
                self.X = numpy.array(lp_problem.solution.get_values())
                self.X[self.enzyme_cols] /= scaling_factor
                self.X[self.process_cols] /= scaling_factor
                self.lambda_ = lp_problem.solution.get_dual_values()
                self.mu_opt = mu_min
            elif exit_flag == INFEASIBLE \
                 or (exit_flag == OPTIMAL_UNSCALED_INFEASIBILITIES):
                mu_max = mu_test
            else:
                print('At mu = ' + str(mu_test) + ': Unknown exit flag '
                      + str(exit_flag) + ' corresponding to status '
                      + lp_problem.solution.get_status_string() + '. '
                      'Interrupting execution')
                return
            # next mu to be tested
            mu_test = (mu_min+mu_max)/2
        print(mu_min)
        return self

    def _setup_lp(self, scaling_factor):
        lp_problem = cplex.Cplex()
    
        # set parameters
        lp_problem.objective.set_sense(lp_problem.objective.sense.minimize)
        lp_problem.parameters.feasopt.tolerance.set(1e-9);
        lp_problem.parameters.simplex.tolerances.feasibility.set(1e-9);
        lp_problem.parameters.simplex.tolerances.optimality.set(1e-9);
        lp_problem.parameters.simplex.tolerances.markowitz.set(0.1);
        lp_problem.parameters.barrier.convergetol.set(1e-9);
        # Threads: the default (0) means that Cplex decides automatically
        # how many threads to use
        # lp_problem.parameters.threads.set(0)
        lp_problem.set_results_stream(None)

        # define columns
        lp_problem.variables.add(obj = self.f, ub = self.UB, lb = self.LB,
                                 names = self.col_names)

        # rescale concentration columns
        scaling = numpy.ones(len(self.col_names))
        scaling[numpy.concatenate([self.enzyme_cols, self.process_cols])] \
            = 1.0/scaling_factor
        S = diags(scaling)
        my_A = self.A * S
        my_Aeq = self.Aeq * S
        
        # prepare inequality and equality constrainsts
        my_A = my_A.tolil()
        my_Aeq = my_Aeq.tolil()
        rows = []
        for nz_ind, data in zip(my_A.rows, my_A.data):
            rows.append(cplex.SparsePair(nz_ind, data))
        for nz_ind, data in zip(my_Aeq.rows, my_Aeq.data):
            rows.append(cplex.SparsePair(nz_ind, data))
        senses = ['L'] * len(self.b) + ['E'] * len(self.beq)
        names = self.ineq_row_names + self.eq_row_names
        rhs = numpy.append(self.b, self.beq)

        lp_problem.linear_constraints.add(lin_expr = rows, rhs = rhs,
                                          senses = senses, names = names)
        return lp_problem

