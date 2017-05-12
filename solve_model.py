import sys

import rba
from rba.post_process import saturating_constraints

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
        solver = model.solve(medium)
        #saturating_constraints(solver)
