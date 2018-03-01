
# python 2/3 compatibility
from __future__ import division, print_function, absolute_import

# global imports

# local imports
from rba.prerba.default_data import DefaultData


class Enzyme(object):
    def __init__(self, reaction, is_transporter):
        self.reaction = reaction
        self.composition = []
        self._is_transporter = is_transporter
        self._initialize_efficiencies()

    def _initialize_efficiencies(self):
        def_activities = DefaultData().activity
        if self._is_transporter:
            self.forward = def_activities.transport_aggregate_id(self.reaction)
            self.backward = def_activities.transport_id
        else:
            self.forward = self.backward = def_activities.efficiency_id
