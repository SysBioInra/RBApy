
from target_vector import TargetVector

class Density(object):
    def __init__(self, max_densities, functions, known_compartments):
        # extract maximal densities
        self.compartments = [md.compartment for md in max_densities]
        self.compartment_indices = [known_compartments.index(md.compartment) \
                                    for md in max_densities]
        self.maximum = TargetVector(max_densities, functions)
