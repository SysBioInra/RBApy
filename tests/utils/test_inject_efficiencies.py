
from __future__ import absolute_import, division, print_function

import rba
from rba.utils.inject_efficiencies import EfficiencyInjecter


def test_enzyme_with_default_efficiencies():
    model = dummy_model()
    EfficiencyInjecter(model, 'tests/utils/efficiencies.tsv')
    enzyme = model.enzymes.enzymes.get_by_id('test_enzyme')
    assert enzyme.forward_efficiency == 'test_enzyme_forward_efficiency'
    assert enzyme.backward_efficiency == 'test_enzyme_backward_efficiency'
    assert model.parameters.functions.get_by_id('test_enzyme_forward_efficiency').parameters[0].value == '0.1'
    assert model.parameters.functions.get_by_id('test_enzyme_backward_efficiency').parameters[0].value == '0.2'


def dummy_model():
    model = rba.RbaModel()
    for e in dummy_enzymes():
        model.enzymes.enzymes.append(e)
    for p in dummy_parameters():
        model.parameters.functions.append(p)
    return model


def dummy_enzymes():
    return (
        rba.xml.Enzyme('test_enzyme', 'test_reaction',
                       'default_efficiency', 'default_efficiency'),
        rba.xml.Enzyme('test_enzyme', 'test_reaction',
                       'default_transporter_efficiency',
                       'default_transporter_efficiency')
    )


def dummy_parameters():
    return (
        rba.xml.Function('default_efficiency', {'CONSTANT': 1}),
        rba.xml.Function('default_transporter_efficiency', {'CONSTANT': 10})
    )
