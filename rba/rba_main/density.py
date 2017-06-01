
from .target_vector import TargetVector

class Density(object):
    def __init__(self, target_densities,
                 known_functions, known_compartments):
        # extract target densities
        self.compartments = [md.compartment for md in target_densities]
        self.compartment_indices = [known_compartments.index(md.compartment) \
                                    for md in target_densities]
        values = []
        self.signs = []
        for target in target_densities:
            if target.value is not None:
                values.append(target.value)
                self.signs.append('E')
            elif target.upper_bound is not None:
                values.append(target.upper_bound)
                self.signs.append('L')
            else:
                raise UserWarning('Density constraint ' + target.compartment
                                  + ': you must specify a value or an upper '
                                  'bound.')
        self.values = TargetVector(values, known_functions)
