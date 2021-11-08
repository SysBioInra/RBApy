
from __future__ import absolute_import, division, print_function

import pytest

import rba
# from rba.utils.inject_efficiencies import EfficiencyInjecter

@pytest.fixture
def model():
    model = rba.RbaModel()
    for e in enzymes():
        model.enzymes.enzymes.append(e)
    for f in functions():
        model.parameters.functions.append(f)
    for a in aggregates():
        model.parameters.aggregates.append(a)
    EfficiencyInjecter(model, 'tests/utils/efficiencies.tsv')
    return model


def enzymes():
    return (
        rba.xml.Enzyme('default_enzyme', 'test_reaction',
                       'default_efficiency', 'default_efficiency'),
        rba.xml.Enzyme('default_transporter', 'test_reaction',
                       'default_transporter_efficiency',
                       'default_transporter_efficiency'),
        rba.xml.Enzyme('parameterized_enzyme', 'test_reaction',
                       'forward_efficiency', 'backward_efficiency'),
        rba.xml.Enzyme('parameterized_transporter', 'test_reaction',
                       'forward_transport', 'backward_transport')
    )


def functions():
    return (
        create_constant('default_efficiency', 1),
        create_constant('default_transporter_efficiency', 2),
        create_constant('forward_efficiency', 3),
        create_constant('backward_efficiency', 4),
        create_constant('transport_factor', 5)
    )


def create_constant(id_, value):
    return rba.xml.Function(id_, 'constant', {'CONSTANT': value})


def aggregates():
    return (
        create_transport_aggregate('forward_transport'),
        create_transport_aggregate('backward_transport')
    )


def create_transport_aggregate(id_):
    result = rba.xml.Aggregate(id_, 'multiplication')
    result.function_references.append(
        rba.xml.FunctionReference('default_transporter_efficiency')
    )
    result.function_references.append(
        rba.xml.FunctionReference('transport_factor')
    )
    return result


@pytest.mark.skip(reason='Tests are outdated')
def test_default_enzyme(model):
    enzyme = model.enzymes.enzymes.get_by_id('default_enzyme')
    assert_standard_efficiency_ids(enzyme)
    assert_efficiencies(model.parameters.functions, enzyme, 0.1, 0.2)


def assert_standard_efficiency_ids(enzyme):
    assert enzyme.forward_efficiency == '{}_forward_efficiency'.format(enzyme.id)
    assert enzyme.backward_efficiency == '{}_backward_efficiency'.format(enzyme.id)


def assert_efficiencies(fns, enzyme, forward, backward):
    assert_constant_parameter(fns, enzyme.forward_efficiency, forward)
    assert_constant_parameter(fns, enzyme.backward_efficiency, backward)


def assert_constant_parameter(functions, id_, expected):
    param = functions.get_by_id(id_).parameters[0]
    assert param.id == 'CONSTANT'
    assert float(param.value) == expected


@pytest.mark.skip(reason='Tests are outdated')
def test_default_transporter(model):
    enzyme = model.enzymes.enzymes.get_by_id('default_transporter')
    assert_standard_efficiency_ids(enzyme)
    assert_efficiencies(model.parameters.functions, enzyme, 0.3, 0.4)


@pytest.mark.skip(reason='Tests are outdated')
def test_parameterized_enzyme(model):
    enzyme = model.enzymes.enzymes.get_by_id('parameterized_enzyme')
    assert enzyme.forward_efficiency == 'forward_efficiency'
    assert enzyme.backward_efficiency == 'backward_efficiency'
    assert_efficiencies(model.parameters.functions, enzyme, 0.5, 0.6)


@pytest.mark.skip(reason='Tests are outdated')
def test_parameterized_transporter(model):
    enzyme = model.enzymes.enzymes.get_by_id('parameterized_transporter')
    assert enzyme.forward_efficiency == 'forward_transport'
    assert enzyme.backward_efficiency == 'backward_transport'
    assert_transport_aggregate(model.parameters, enzyme, 0.7, 0.8)


def assert_transport_aggregate(params, enzyme, forward, backward):
    assert_aggregate_efficiency(params, enzyme.forward_efficiency,
                                '{}_forward_efficiency'.format(enzyme.id),
                                forward)
    assert_aggregate_efficiency(params, enzyme.backward_efficiency,
                                '{}_backward_efficiency'.format(enzyme.id),
                                backward)


def assert_aggregate_efficiency(params, agg_id, fn_id, expected):
    agg = params.aggregates.get_by_id(agg_id)
    assert len(agg.function_references) == 2
    assert agg.function_references[0].function == fn_id
    assert agg.function_references[1].function == 'transport_factor'
    assert_constant_parameter(params.functions, fn_id, expected)
