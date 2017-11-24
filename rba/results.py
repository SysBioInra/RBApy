

from __future__ import absolute_import, division, print_function

import itertools


class Results(object):

    def __init__(self, variables, dual_values, model):
        self.variables = variables
        self.dual_values = dual_values
        self.reactions = {r.id: reaction_string(r)
                          for r in model.metabolism.reactions}
        self.enzymes = [e.id for e in model.enzymes.enzymes]
        self.processes = [p.id for p in model.processes.processes]
        boundary_metabolites = set(m.id for m in model.metabolism.species
                                   if m.boundary_condition)
        self.boundary_reactions = []
        for r in model.metabolism.reactions:
            if any(m.species in boundary_metabolites
                   for m in itertools.chain(r.reactants, r.products)):
                self.boundary_reactions.append(r.id)

    def reaction_fluxes(self):
        return {i: self.variables[i] for i in self.reactions}

    def enzyme_concentrations(self):
        return {i: self.variables.get(i + '_enzyme', 0) for i in self.enzymes}

    def process_machinery_concentrations(self):
        return {i: self.variables[i + '_machinery'] for i in self.processes}

    def sorted_boundary_fluxes(self):
        """
        Return all nonzero import reactions.

        Parameters
        ----------
        solver : RbaSolver
            Solver storing RBA problem solution.

        Returns
        -------
        out : List of (reaction_id, flux) tuples
            List is sorted by reactions with higher flux first.

        """
        return self.sorted_fluxes(self.boundary_reactions)

    def sorted_fluxes(self, ids=None):
        if not ids:
            ids = list(self.reactions.keys())
        result = []
        for reaction in ids:
            flux = self.variables[reaction]
            if flux != 0:
                result.append((self.reactions[reaction], flux))
        result.sort(key=lambda x: abs(x[1]), reverse=True)
        return result

    # def saturating_constraints(self):
    #     """
    #     Print relevant saturated constraints.
    #
    #     Check for (i) saturated reaction lower bounds and upper bounds,
    #     (ii) saturated density constraints.
    #
    #     Parameters
    #     ----------
    #     solver : RbaSolver
    #         Solver used to solve RBA problem.
    #
    #     """
    #     # find saturating lb or ub
    #     for c in solver.reaction_cols:
    #         nu = solver.X[c]
    #         ub = solver.UB[c]
    #         lb = solver.LB[c]
    #         if lb != 0 and abs((nu-lb)/lb) < 0.1:
    #             print('{}: {} vs LB {}'.format(solver.col_names[c], nu, lb))
    #         if ub != 0 and abs((nu-ub)/ub) < 0.1:
    #             print('{}: {} vs LB {}'.format(solver.col_names[c], nu, ub))
    #
    #     # find saturating density constraints
    #     density_rows = [i for i, r_name in enumerate(solver.row_names)
    #                     if 'density' in r_name]
    #     b_opt = solver.A.dot(solver.X)
    #     for r in density_rows:
    #         # only display info if constraint is saturated > 95%
    #         if solver.b[r] == 0 or (b_opt[r] / solver.b[r] < 0.95):
    #             continue
    #         # find critical enzymes along this constraint
    #         keep = 10
    #         print('\n' + solver.row_names[r]
    #               + ' (%g vs %g)' % (solver.b[r], b_opt[r])
    #               + '. Top {}:'.format(keep))
    #         enzyme_density = (solver.A.toarray()[r, :] * solver.X) / solver.b[r]
    #         order = numpy.argsort(enzyme_density)[::-1][:keep]
    #         for i in order:
    #             print('{} {}%'.format(solver.col_names[i], 100*enzyme_density[i]))
    #
    #     # find non saturated enzymes
    #     output = []
    #     for c in solver.enzyme_cols:
    #         # find rows containing forward and backward capacity
    #         enzyme = solver.col_names[c]
    #         f = solver.row_names.index(enzyme + '_forward_capacity')
    #         b = solver.row_names.index(enzyme + '_backward_capacity')
    #         if (b_opt[f] < -1e-10) and (b_opt[b] < -1e-10):
    #             output.append(enzyme + ': (%g, %g)' % (b_opt[f], b_opt[b]))
    #     if len(output) > 0:
    #         print('\nNon saturated enzymes (forward, backward):\n'
    #               + '\n'.join(output))


def reaction_string(reaction):
    reactants = ' + '.join(['{} {}'.format(r.stoichiometry, r.species)
                            for r in reaction.reactants])
    products = ' + '.join(['{} {}'.format(p.stoichiometry, p.species)
                           for p in reaction.products])
    return reaction.id + ': ' + reactants + ' -> ' + products
