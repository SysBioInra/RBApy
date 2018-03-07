
# python 2/3 compatibility
from __future__ import division, print_function, absolute_import

# global imports

# local imports
from rba.prerba.default_data import DefaultActivity


class Enzyme(object):
    DEF_ACTIVITIES = DefaultActivity()

    def __init__(self, reaction, is_transporter):
        self.reaction = reaction
        self.gene_association = []
        self.composition = []
        self.imported_metabolites = []
        self.is_transporter = is_transporter
        self.initialize_efficiencies()

    def initialize_efficiencies(self):
        if self.is_transporter:
            if self.imported_metabolites:
                self.forward = self.DEF_ACTIVITIES.transport_aggregate_id(
                    self.reaction
                )
            else:
                self.forward = self.DEF_ACTIVITIES.transport_id
            self.backward = self.DEF_ACTIVITIES.transport_id
        else:
            self.forward = self.backward = self.DEF_ACTIVITIES.efficiency_id
