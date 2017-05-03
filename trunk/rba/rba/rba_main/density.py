
from target_vector import TargetVector

class Density(object):
    def __init__(self, data, functions):
        # extract maximal densities
        max_densities = data.parameters.maximal_densities
        known_compartments = [c.id for c in data.metabolism.compartments]
        self.compartments = [md.compartment for md in max_densities]
        self.compartment_indices = [known_compartments.index(md.compartment) \
                                    for md in max_densities]
        self.maximum = TargetVector(max_densities, functions)
