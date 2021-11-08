#!/usr/bin/env python3
"""Build RBA XML files from SBML."""

import argparse
import sys

import rba


def main():
    parser = argparse.ArgumentParser(description='Generate an RNA model')
    parser.add_argument('parameter_file', metavar='parameter-file', type=str,
                        help='Path to parameter file')
    parser.add_argument('--model-dir', type=str, default=None,
                        help='Directory to save model (default: parent directory of the the parameter file)')
    parser.add_argument('--verbose', action='store_true',
                        help='Directory to save the results; defaut: the directory for the model.')

    args = parser.parse_args()

    if args.verbose:
        print('Building model ...')
    try:
        model = rba.RbaModel.from_data(args.parameter_file, verbose=args.verbose)
        model.write(args.model_dir)
    except UserWarning as error:
        raise SystemExit(str(error))
    if args.verbose:
        print('done')


if __name__ == '__main__':
    main()
