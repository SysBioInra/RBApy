"""Module providing simple post processing functions."""

import numpy


def nonzero_imports(solver):
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
    S = solver._blocks.S.tocsc()
    import_indices = [i for i in range(S.shape[1]) if S[:, i].nnz == 1]
    import_reactions = [solver.col_names[i] for i in import_indices]
    import_fluxes = solver.X[import_indices]
    nz = numpy.nonzero(import_fluxes)[0]
    nz_import_reactions = [import_reactions[i] for i in nz]
    imports = list(zip(nz_import_reactions, import_fluxes[nz]))
    imports.sort(key=lambda x: x[1], reverse=True)
    return imports


def saturating_constraints(solver):
    """
    Print relevant saturated constraints.

    Check for (i) saturated reaction lower bounds and upper bounds,
    (ii) saturated density constraints.

    Parameters
    ----------
    solver : RbaSolver
        Solver used to solve RBA problem.

    """
    # find saturating lb or ub
    for c in solver.reaction_cols:
        nu = solver.X[c]
        ub = solver.UB[c]
        lb = solver.LB[c]
        if lb != 0 and abs((nu-lb)/lb) < 0.1:
            print('{}: {} vs LB {}'.format(solver.col_names[c], nu, lb))
        if ub != 0 and abs((nu-ub)/ub) < 0.1:
            print('{}: {} vs LB {}'.format(solver.col_names[c], nu, ub))

    # find saturating density constraints
    density_rows = [i for i, r_name in enumerate(solver.row_names)
                    if 'density' in r_name]
    b_opt = solver.A.dot(solver.X)
    for r in density_rows:
        # only display info if constraint is saturated > 95%
        if solver.b[r] == 0 or (b_opt[r] / solver.b[r] < 0.95):
            continue
        # find critical enzymes along this constraint
        keep = 10
        print('\n' + solver.row_names[r]
              + ' (%g vs %g)' % (solver.b[r], b_opt[r])
              + '. Top {}:'.format(keep))
        enzyme_density = (solver.A.toarray()[r, :] * solver.X) / solver.b[r]
        order = numpy.argsort(enzyme_density)[::-1][:keep]
        for i in order:
            print('{} {}%'.format(solver.col_names[i], 100*enzyme_density[i]))

    # find non saturated enzymes
    output = []
    for c in solver.enzyme_cols:
        # find rows containing forward and backward capacity
        enzyme = solver.col_names[c]
        f = solver.row_names.index(enzyme + '_forward_capacity')
        b = solver.row_names.index(enzyme + '_backward_capacity')
        if (b_opt[f] < -1e-10) and (b_opt[b] < -1e-10):
            output.append(enzyme + ': (%g, %g)' % (b_opt[f], b_opt[b]))
    if len(output) > 0:
        print('\nNon saturated enzymes (forward, backward):\n'
              + '\n'.join(output))
