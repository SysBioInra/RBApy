"""
Inject efficiencies into existing model.

New efficiencies are read from a tab-separated CSV file.
Each line of the file must have one of the 4 following formats
(all spaces are really tabs):
 - default_efficiency <default_efficiency_value>
 - default_transporter_efficiency <default_transporter_efficiency_value>
 - <enzyme_id> <forward_efficiency_constant> <backward_efficiency_constant>
 - <enzyme_id> <forward/backward> <fn_type> [<param_name> <param_value>]^n

The first two formats enable to update default efficiencies.
The third format updates one enzyme's forward and backward efficiencies
with constant parameters.
The last format updates one enzyme's forward or backward efficiency with
an arbitrary parameter function.
Note that for transporters, only the base efficiency is modified, not the
Michaelis-Menten term!
More specifically, enzyme-specific updates occur as follows:
 - current efficiency is a default efficiency: it is replaced by a new
parameter function with id '<enzyme_id>_<forward/backward>_efficiency'.
 - current efficiency is a parameter function: we keep the current
parameter id but replace it with the new efficiency.
 - current efficiency is a parameter aggregate (e.g. enzyme is a transporter):
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


def set_efficiencies(rba_model, efficiency_file):
    EfficiencySetter(rba_model, efficiency_file)


class EfficiencySetter(object):
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
        if efficiency.forward:
            enzyme.forward_efficiency = self._parameter_updater.update(
                enzyme.forward_efficiency, efficiency.forward
            ).id
        if efficiency.backward:
            enzyme.backward_efficiency = self._parameter_updater.update(
                enzyme.backward_efficiency, efficiency.backward
            ).id


class ParameterUpdater(object):
    def __init__(self, parameters):
        self._parameters = parameters

    def update(self, param_id, efficiency_fn):
        self._efficiency_fn = efficiency_fn
        new_fn = self._update_function(param_id)
        if new_fn:
            return new_fn
        return self._update_aggregate(param_id)

    def _update_function(self, fn_id):
        if not self._is_default(fn_id):
            current_fn = self._parameters.functions.get_by_id(fn_id)
            if not current_fn:
                return None
            self._efficiency_fn.id = current_fn.id
        self._parameters.functions.append(self._efficiency_fn)
        return self._efficiency_fn

    def _is_default(self, parameter_id):
        return parameter_id in ['default_efficiency',
                                'default_transporter_efficiency']

    def _update_aggregate(self, agg_id):
        result = self._parameters.aggregates.get_by_id(agg_id)
        modifiable_ref = self._find_modifiable_function_reference(result)
        if modifiable_ref:
            new_fn = self._update_function(modifiable_ref.function)
            if new_fn:
                modifiable_ref.function = new_fn.id
                return result
        raise UserWarning('Parameter {} could not be updated'.format(agg_id))

    def _find_modifiable_function_reference(self, aggregate):
        for fn_ref in aggregate.function_references:
            if (fn_ref.function == self._efficiency_fn.id
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
                tokens = line.strip().split('\t')
                if (not self._parse_default(tokens)
                        and not self._parse_constants(tokens)
                        and not self._parse_function(tokens)):
                    raise UserWarning('Invalid line: ' + line)

    def _parse_default(self, tokens):
        try:
            id_, value = tokens
        except ValueError:
            return False
        if id_ == 'default_efficiency':
            self.default = value
        elif id_ == 'default_transporter_efficiency':
            self.transporter_default = value
        else:
            return False
        return True

    def _parse_constants(self, tokens):
        try:
            id_, forward, backward = tokens
        except ValueError:
            return False
        eff = self._efficiencies.setdefault(id_, EfficiencyParameters(id_))
        eff.set_constants(forward, backward)
        return True

    def _parse_function(self, tokens):
        try:
            id_, sense, fn_type = tokens[0:3]
            parameters = tokens[3:]
        except ValueError:
            return False
        if (len(parameters) == 0) or (len(parameters) % 2) != 0:
            return False
        parameters = dict(zip(parameters[::2], parameters[1::2]))
        eff = self._efficiencies.setdefault(id_, EfficiencyParameters(id_))
        if sense == 'forward':
            eff.set_forward(fn_type, parameters)
        elif sense == 'backward':
            eff.set_backward(fn_type, parameters)
        else:
            return False
        return True

    def __getitem__(self, id_):
        return self._efficiencies.get(id_)


class EfficiencyParameters(object):
    def __init__(self, enzyme_id):
        self.id = enzyme_id
        self.forward = None
        self.backward = None

    def set_constants(self, forward, backward):
        self.forward = self._create_constant(
            self._standard_id('forward'), forward
        )
        self.backward = self._create_constant(
            self._standard_id('backward'), backward
        )

    def _create_constant(self, id_, efficiency):
        return rba.xml.Function(id_, 'constant', {'CONSTANT': efficiency})

    def _standard_id(self, sense):
        return '{}_{}_efficiency'.format(self.id, sense)

    def set_forward(self, fn_type, parameters):
        self.forward = rba.xml.Function(
            self._standard_id('forward'), fn_type, parameters
        )

    def set_backward(self, fn_type, parameters):
        self.backward = rba.xml.Function(
            self._standard_id('backward'), fn_type, parameters
        )
