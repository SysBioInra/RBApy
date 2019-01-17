"""Solve standard RBA problem."""

# python 2/3 compatibility
from __future__ import division, print_function

# global imports
import sys
import os.path

# package imports
import rba


def main():
    if len(sys.argv) < 2:
        print('Please provide path to directory containing xml files.')
    else:
        xml_dir = sys.argv[1]
        if len(sys.argv) >= 3:
            output_dir = sys.argv[2]
        else:
            output_dir = xml_dir

        # load model, build matrices and solve
        model = rba.RbaModel.from_xml(xml_dir)
        results = model.solve()
	
        print('Optimal growth rate is {}.'.format(results.mu_opt))
        results.write(output_dir)
	#results.write_fluxes(xml_dir+'fluxEcoli',
        #           file_type='json',
        #           merge_isozyme_reactions=True,
        #           only_nonzero=True,
        #           remove_prefix=True)

	
        # results.print_main_transport_reactions()


if __name__ == '__main__':
    main()
