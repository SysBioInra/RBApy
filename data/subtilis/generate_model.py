"""Module generating subtilis model."""

# python 2 compatibility
from __future__ import absolute_import, division, print_function

# imports
import sys
from os import path
import math

sys.path = [path.join(sys.path[0], '../..')] + sys.path
import rba  # noqa

old_data = 'data/subtilis_ref/old_data/'


def main():
    builder = rba.ModelBuilder('data/subtilis/params.in')
    subtilis = builder.build_model()
    subtilis.medium = reference_medium(subtilis)
    add_enzymatic_activities(subtilis.enzymes, subtilis.parameters, 'medium_2')
    apply_old_stoichiometries(subtilis.enzymes)
    add_flagella_constraint(subtilis)
    subtilis.write()


def add_flagella_constraint(subtilis):
    subtilis.targets.target_groups.append(flagella_activation())
    for fn in flagella_activation_functions():
        subtilis.parameters.functions.append(fn)
    subtilis.parameters.aggregates.append(flagella_activation_aggregate())


def flagella_activation():
    target_group = rba.xml.TargetGroup('flagella_activation')
    target = rba.xml.TargetReaction('Th')
    target.value = 'flagella_proton_flux'
    target_group.reaction_fluxes.append(target)
    return target_group


def flagella_activation_functions():
    return [
        rba.xml.Function('flagella_speed', 'constant', {'CONSTANT': 5.81}),
        rba.xml.Function('flagella_h_consumption', 'constant',
                         {'CONSTANT': 0.9415}),
        rba.xml.Function('number_flagella', 'linear',
                         {'LINEAR_COEF': 4.5197, 'LINEAR_CONSTANT': 3.7991,
                          'X_MIN': 0.25, 'X_MAX': 1.6,
                          'Y_MIN': float('-Inf'), 'Y_MAX': float('Inf')})
        ]


def flagella_activation_aggregate():
    aggregate = rba.xml.Aggregate('flagella_proton_flux', 'multiplication')
    aggregate.function_references.append(
        rba.xml.FunctionReference('number_flagella')
    )
    aggregate.function_references.append(
        rba.xml.FunctionReference('flagella_speed')
    )
    aggregate.function_references.append(
        rba.xml.FunctionReference('flagella_h_consumption')
    )
    return aggregate


def add_enzymatic_activities(enzymes, parameters, medium):
    with open(old_data + 'catalytic_activity.csv', 'r') as f:
        activity = read_activities(f, medium)
    for enzyme in enzymes.enzymes:
        name = enzyme.reaction + '_base_efficiency'
        parameters.functions.append(
            rba.xml.Function(name, 'constant',
                             activity[old_name(enzyme.reaction)])
        )
        if is_transport_function(enzyme.forward_efficiency):
            transport_agg = parameters.aggregates.get_by_id(
                enzyme.forward_efficiency
            )
            transport_agg.function_references[0].function = name
        else:
            enzyme.forward_efficiency = name
        enzyme.backward_efficiency = name


def is_transport_function(attribute):
    return not(
        attribute in ['default_efficiency', 'default_transporter_efficiency']
    )


def add_zero_cost_flags(enzymes):
    with open(old_data + 'zero_cost.csv', 'r') as f:
        zero_cost = []
        for line in f:
            zero_cost += line.rstrip('\n').split('\t')
    for enzyme in enzymes.enzymes:
        id_ = old_name(enzyme.id)
        if id_ in zero_cost:
            enzyme.zero_cost = True
            zero_cost.remove(id_)


def read_activities(f, medium):
    fn_type = None
    result = {}
    for line in f:
        token = line.rstrip('\n').split('\t')
        if token[0] == 'function':
            # format 'function id type'
            if token[1] == medium:
                fn_type = token[2]
        else:
            # format 'reaction_id fn_id param1_id param1_value ... paramn_value
            if token[1] == medium:
                reaction = token[0]
                params = iter(token[2:])
                result[reaction] = {
                    id_: value for id_, value in zip(params, params)
                    }
    if not fn_type:
        raise UserWarning('Could not retrieve medium {}.'.format(medium))
    return result


def old_name(new_name):
    if new_name == 'R_maintenance_atp':
        return 'Eatpm'
    else:
        return new_name.rsplit('_', 1)[0]


def apply_old_stoichiometries(enzymes):
    with open(old_data + 'stoichiometry.csv', 'r') as f:
        data = {}
        for line in f:
            [species, sto] = line.rstrip('\n').split('\t')
            data[species] = float(sto)
    for enzyme in enzymes.enzymes:
        for sr in enzyme.machinery_composition.reactants:
            try:
                sr.stoichiometry = data[sr.species]
            except KeyError:
                pass


def reference_medium(subtilis):
    with open(old_data + 'medium.csv', 'r') as f:
        old_medium = dict.fromkeys(subtilis.medium, 0)
        for line in f:
            met, conc = line.rstrip('\n').split('\t')
            old_medium[met[:-2]] = conc
    return old_medium


if __name__ == "__main__":
    main()
