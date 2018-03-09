"""
Inject efficiencies into existing model.

We parse the CSV file with new enzyme efficiencies. There are 3 cases:
 - current efficiency is a default efficiency: a new parameter
function with id '<enzyme_id>_<forward/backward>_efficiency' is created.
 - current efficiency points to a parameter function: we keep the current
parameter id but replace it with the new efficiency.
 - current efficiency points to a parameter aggregate: we assume that the
 enzyme is a transporter.

If the enzyme is assumed to be a transporter, there are again 3 cases:
 - one of the function of the aggregate points to a default efficiency:
 we create a new parameter function with id
 '<enzyme_id>_<forward/backward>_efficiency' that replaces
 the default efficiency in the aggregate.
 - one of the function of the aggregate points to
 '<enzyme_id>_<forward/backward>_efficiency': we replace this function with
 the new efficiency.
 - none of the above apply: we print a warning and do not modify anything.
"""

from __future__ import absolute_import, division, print_function

import rba


class EfficiencyInjecter(object):
    def __init__(self, rba_model, efficiency_file):
        self.model = rba_model
        self._parameter_updater = ParameterUpdater(self.model.parameters)
        self._efficiencies = EfficiencyFile(efficiency_file)
        self._update_default_efficiencies()
        self._inject_efficiencies()

    def _update_default_efficiencies(self):
        if self._efficiencies.default:
            self.model.parameters.functions.get_by_id(
                'default_efficiency'
            ).value = self._efficiencies.default
        if self._efficiencies.transporter_default:
            self.model.parameters.functions.get_by_id(
                'default_transporter_efficiency'
            ).value = self._efficiencies.transporter_default

    def _inject_efficiencies(self):
        for enzyme in self.model.enzymes.enzymes:
            new_efficiency = self._efficiencies[enzyme.id]
            if new_efficiency:
                self._change_efficiency(enzyme, new_efficiency)

    def _change_efficiency(self, enzyme, efficiency):
        enzyme.forward_efficiency = self._parameter_updater.update(
            enzyme.forward_efficiency,
            self._standard_id(enzyme.id, 'forward'),
            efficiency[0]
        ).id
        enzyme.backward_efficiency = self._parameter_updater.update(
            enzyme.backward_efficiency,
            self._standard_id(enzyme.id, 'backward'),
            efficiency[1]
        ).id

    def _standard_id(self, enzyme_id, sense):
        return '{}_{}_efficiency'.format(enzyme_id, sense)


class ParameterUpdater(object):
    def __init__(self, parameters):
        self._parameters = parameters

    def update(self, param_id, standard_id, efficiency):
        self._standard_id = standard_id
        self._efficiency = efficiency
        new_fn = self._update_function(param_id)
        if new_fn:
            return new_fn
        return self._update_aggregate(param_id)

    def _update_function(self, fn_id):
        if self._is_default(fn_id):
            new_fn = self._create_function(self._standard_id)
            self._parameters.functions.append(new_fn)
            return new_fn
        current_fn = self._parameters.functions.get_by_id(fn_id)
        if current_fn:
            current_fn = self._create_function(current_fn.id)
            return current_fn
        return None

    def _is_default(self, parameter_id):
        return parameter_id in ['default_efficiency',
                                'default_transporter_efficiency']

    def _create_function(self, id_):
        return rba.xml.Function(id_, 'constant',
                                {'CONSTANT': self._efficiency})

    def _update_aggregate(self, agg_id):
        result = self._parameters.aggregates.get_by_id(agg_id)
        modifiable_ref = self._find_modifiable_function_reference(result)
        if modifiable_ref:
            new_fn = self._update_function(modifiable_ref.function)
            if new_fn:
                modifiable_ref.function = new_fn.id
                return result
        raise UserWarning('Aggregate could not be updated')

    def _find_modifiable_function_reference(self, aggregate):
        for fn_ref in aggregate.function_references:
            if (fn_ref.function == self._standard_id
                    or self._is_default(fn_ref.function)):
                return fn_ref
        return None


class EfficiencyFile(object):
    def __init__(self, filename):
        self.default = None
        self.transporter_default = None
        self._efficiencies = {}
        with open(filename, 'r') as input_stream:
            for line in input_stream:
                try:
                    self._parse_efficiency(line)
                except ValueError:
                    self._parse_default(line)

    def _parse_efficiency(self, line):
        id_, forward, backward = line.strip().split('\t')
        self._efficiencies[id_] = (forward, backward)

    def _parse_default(self, line):
        id_, value = line.strip().split('\t')
        if id_ == 'default_efficiency':
            self.default = value
        elif id_ == 'default_transporter_efficiency':
            self.transporter_default = value
        else:
            raise UserWarning('Invalid line: ' + line)

    def __getitem__(self, id_):
        return self._efficiencies.get(id_)
