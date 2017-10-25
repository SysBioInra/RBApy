"""Solve standard RBA problem."""

# python 2/3 compatibility
from __future__ import division, print_function

# global imports
import sys

# package imports
import rba


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('Please provide path to directory containing xml files.')
    else:
        xml_dir = sys.argv[1]
        if len(sys.argv) >= 3:
            medium = sys.argv[2]
        else:
            medium = 'default'

        model = rba.RbaModel.from_xml(xml_dir)
        results = model.solve(medium)
