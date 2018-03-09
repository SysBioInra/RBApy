
from __future__ import absolute_import, division, print_function

import pytest

import rba
from rba.utils.inject_efficiencies import EfficiencyInjecter


@pytest.fixture
def model():
    model = rba.RbaModel()
    for e in enzymes():
        model.enzymes.enzymes.append(e)
    for f in functions():
        model.parameters.functions.append(f)
    EfficiencyInjecter(model, 'tests/utils/efficiencies.tsv')
    return model


def enzymes():
    return (
        rba.xml.Enzyme('default_enzyme', 'test_reaction',
                       'default_efficiency', 'default_efficiency'),
        rba.xml.Enzyme('default_transporter', 'test_reaction',
                       'default_transporter_efficiency',
                       'default_transporter_efficiency'),
        rba.xml.Enzyme('enzyme_parameter_efficiencies', 'test_reaction',
                       'forward_efficiency', 'backward_efficiency')
    )


def functions():
    return (
        create_constant('default_efficiency', 1),
        create_constant('default_transporter_efficiency', 2),
        create_constant('forward_efficiency', 3),
        create_constant('backward_efficiency', 4)
    )


def create_constant(id_, value):
    return rba.xml.Function(id_, 'constant', {'CONSTANT': value})


def test_default_enzyme(model):
    enzyme = model.enzymes.enzymes.get_by_id('default_enzyme')
    assert_standard_efficiency_ids(enzyme)
    assert_efficiencies(model.parameters.functions, enzyme, 0.1, 0.2)


def assert_standard_efficiency_ids(enzyme):
    assert enzyme.forward_efficiency == '{}_forward_efficiency'.format(enzyme.id)
    assert enzyme.backward_efficiency == '{}_backward_efficiency'.format(enzyme.id)


def assert_efficiencies(fns, enzyme, forward, backward):
    assert get_constant_parameter(fns, enzyme.forward_efficiency) == forward
    assert get_constant_parameter(fns, enzyme.backward_efficiency) == backward


def get_constant_parameter(functions, id_):
    param = functions.get_by_id(id_).parameters[0]
    assert param.id == 'CONSTANT'
    return float(param.value)


def test_default_transporter(model):
    enzyme = model.enzymes.enzymes.get_by_id('default_transporter')
    assert_standard_efficiency_ids(enzyme)
    assert_efficiencies(model.parameters.functions, enzyme, 0.3, 0.4)


def test_enzyme_parameter_efficiencies(model):
    enzyme = model.enzymes.enzymes.get_by_id('enzyme_parameter_efficiencies')
    assert enzyme.forward_efficiency == 'forward_efficiency'
    assert enzyme.backward_efficiency == 'backward_efficiency'
    assert_efficiencies(model.parameters.functions, enzyme, 0.5, 0.6)
