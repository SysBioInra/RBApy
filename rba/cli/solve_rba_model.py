#!/usr/bin/env python3
"""Solve standard RBA problem."""

# python 2/3 compatibility
from __future__ import division, print_function

# global imports
import argparse
import sys

# package imports
import rba


def main():
    parser = argparse.ArgumentParser(description='Solve an RNA model')
    parser.add_argument('model_dir', metavar='model-dir', type=str,
                        help='Directory for model')
    parser.add_argument('--lp-solver', type=str, default=None,
                        help=(
                            'LP solver (`cplex`, `glpk`, `swiglpk`, ``gurobi``, or `scipy`; '
                            'default: `cplex` if installed, otherwise `gurobi` if installed, otherwise `glpk`)'
                        ))
    parser.add_argument('--mu-min', type=float, default=0.,
                        help='Minimum μ to check; default: 0.0.')
    parser.add_argument('--mu-max', type=float, default=2.5,
                        help='Maximum μ to check; default: 2.5.')
    parser.add_argument('--bissection-tol', type=float, default=1e-4,
                        help='Tolerance for bissection; default: 1e-4.')
    parser.add_argument('--max-bissection-iters', type=int, default=None,
                        help='Maximum number of iterations for bissection; default: None.')
    parser.add_argument('--output-dir', type=str, default=None,
                        help='Directory to save the results; defaut: the directory for the model.')
    parser.add_argument('--verbose', action='store_true',
                        help='Directory to save the results; defaut: the directory for the model.')

    args = parser.parse_args()

    # load model, build matrices and solve
    if args.verbose:
        print('Reading model from XML files ...', end='')
        sys.stdout.flush()
    model = rba.RbaModel.from_xml(args.model_dir)
    if args.verbose:
        print(' done')

    if args.verbose:
        print('Solving RBA model ...')
    results = model.solve(lp_solver=args.lp_solver,
                          mu_min=args.mu_min,
                          mu_max=args.mu_max,
                          bissection_tol=args.bissection_tol,
                          max_bissection_iters=args.max_bissection_iters,
                          verbose=args.verbose)
    if args.verbose:
        print('done')

    print('Optimal growth rate is {}.'.format(results.mu_opt))
    results.write(args.output_dir or args.model_dir)


if __name__ == '__main__':
    main()
