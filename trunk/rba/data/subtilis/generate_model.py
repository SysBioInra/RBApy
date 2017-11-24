"""Module generating subtilis model."""

# python 2 compatibility
from __future__ import absolute_import, division, print_function

# imports
import os.path
import sys
import math

sys.path.append(os.path.join(sys.path[0], '../..'))
import rba  # noqa

old_data = 'data/subtilis_ref/old_data/'


def flagella_activation():
    process = rba.xml.Process('P_FLAGELLA', 'Flagella activation')
    target = rba.xml.TargetReaction('Th')
    target.value = 'flagella_proton_flux'
    process.targets.reaction_fluxes.append(target)
    return process


def flagella_activation_functions():
    return [
        Function('flagella_speed', 'constant', {'CONSTANT': 5.81}),
        Function('flagella_h_consumption', 'constant',
                 {'CONSTANT': 0.9415}),
        Function('number_flagella', 'linear',
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
        reaction = old_name(enzyme.enzymatic_activity.reaction)
        pass

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
    params = {}
    for line in f:
        token = line.rstrip('\n').split('\t')
        if token[0] == 'function':
            # format 'function id type'
            if token[1] == medium:
                fn_type = token[2]
        else:
            # format 'reaction_id fn_id param1_id param1_value ... paramn_value
            fn = token[1]
            if fn == 'medium':
                reaction = token[0]
                params = iter(token[2:])
                params[reaction] = {
                    id_: value for id_, value in zip(params, params)
                    }
    if not fn_type:
        raise UserWarning('Could not retrieve medium {}.'.format(medium))
    return params


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
                print(sr.species)


if __name__ == "__main__":
    # generate base XML files with pipeline
    subtilis = rba.prerba.PreRba('data/subtilis_test/params.in')

    # add flagella constraint
    subtilis.model.processes.processes.append(flagella_activation())
    for fn in flagella_activation_functions():
        subtilis.model.parameters.functions.append(fn)
    subtilis.model.parameters.aggregates.append(
        flagella_activation_aggregate()
        )

    # add enzymatic activities
    add_enzymatic_activities(subtilis.model.enzymes,
                             subtilis.model.parameters, 'medium_2')

    # add zero_cost flags
    # add_zero_cost_flags(subtilis.model.enzymes)

    # apply old stoichiometries
    apply_old_stoichiometries(subtilis.model.enzymes)

    # set medium to original medium
    with open(old_data + 'medium.csv', 'r') as f:
        old_medium = {}
        for line in f:
            met, conc = line.rstrip('\n').split('\t')
            old_medium[met[:-2]] = conc
    subtilis.model.medium = dict.fromkeys(subtilis.model.medium, 0)
    for met, conc in old_medium.items():
        subtilis.model.medium[met] = conc

    # write xml files
    subtilis.model.write_files()
