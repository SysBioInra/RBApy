"""
Module defining enzyme-specific classes used for RBA XML structures.
"""

# python 2/3 compatiblity
from __future__ import division, print_function

# global imports
from lxml import etree

# local imports
from rba.xml.common import (get_unique_child, is_true, ListOf,
                            Function, MachineryComposition,
                            Parameter, ListOfParameters)

__all__ = ['RbaEnzymes', 'ListOfEfficiencyFunctions', 'Enzyme', 'ListOfEnzymes',
           'EnzymaticActivity', 'TransporterEfficiency',
           'EnzymeEfficiency', 'ListOfEnzymeEfficiencies']


class RbaEnzymes(object):
    def __init__(self):
        self.efficiency_functions = ListOfEfficiencyFunctions()
        self.enzymes = ListOfEnzymes()

    def write(self, output_stream, doc_name = 'RBAEnzymes'):
        root = etree.Element(doc_name)
        root.extend([self.efficiency_functions.to_xml_node(),
                     self.enzymes.to_xml_node()])
        etree.ElementTree(root).write(output_stream, pretty_print=True)

    @classmethod
    def from_file(cls, input_stream):
        node = etree.ElementTree(file=input_stream).getroot()
        result = cls()
        n = get_unique_child(node, ListOfEfficiencyFunctions.tag)
        result.efficiency_functions = ListOfEfficiencyFunctions.from_xml_node(n)
        n = get_unique_child(node, ListOfEnzymes.tag)
        result.enzymes = ListOfEnzymes.from_xml_node(n)
        return result

class ListOfEfficiencyFunctions(ListOf):
    tag = 'listOfEfficiencyFunctions'
    list_element = Function

class Enzyme(object):
    tag = 'enzyme'
    
    def __init__(self, reaction, zero_cost=False):
        if zero_cost is None: zero_cost = False
        self.id = reaction
        self.zero_cost = zero_cost
        self.machinery_composition = MachineryComposition()
        self.enzymatic_activity = EnzymaticActivity(reaction)

    def to_xml_node(self):
        result = etree.Element(self.tag)
        result.set('id', self.id)
        result.set('zeroCost', str(self.zero_cost).lower())
        if not(self.machinery_composition.is_empty()):
            result.extend([self.machinery_composition.to_xml_node()])
        result.extend([self.enzymatic_activity.to_xml_node()])
        return result

    @classmethod
    def from_xml_node(cls, node):
        result = cls(node.get('id'), is_true(node.get('zeroCost')))
        n = get_unique_child(node, MachineryComposition.tag, False)
        if n is not None:
            result.machinery_composition = MachineryComposition.from_xml_node(n)
        n = get_unique_child(node, EnzymaticActivity.tag)
        result.enzymatic_activity = EnzymaticActivity.from_xml_node(n)
        return result

class ListOfEnzymes(ListOf):
    tag = 'listOfEnzymes'
    list_element = Enzyme

class EnzymaticActivity(object):
    tag = 'enzymaticActivity'
    
    def __init__(self, reaction):
        self.reaction = reaction
        self.enzyme_efficiencies = ListOfEnzymeEfficiencies()
        self.transporter_efficiency = TransporterEfficiency()

    def to_xml_node(self):
        result = etree.Element(self.tag)
        result.set('reaction', self.reaction)
        # optional subelements
        opt = [self.enzyme_efficiencies, self.transporter_efficiency]
        result.extend([e.to_xml_node() for e in opt if not(e.is_empty())])
        return result

    @classmethod
    def from_xml_node(cls, node):
        reaction = node.get('reaction')
        result = cls(reaction)
        n = get_unique_child(node, ListOfEnzymeEfficiencies.tag)
        result.enzyme_efficiencies = ListOfEnzymeEfficiencies.from_xml_node(n)
        n = get_unique_child(node, TransporterEfficiency.tag, False)
        if n is not None: result.transporter_efficiency = \
           TransporterEfficiency.from_xml_node(n)
        return result

class EnzymeEfficiency(object):
    tag = 'enzymeEfficiency'
    
    def __init__(self, function, parameters=None):
        self.function = function
        self.parameters = ListOfParameters()
        if parameters:
            for key,value in parameters.items():
                self.parameters.append(Parameter(key, value))

    def to_xml_node(self):
        result = etree.Element(self.tag)
        result.set('function', self.function)
        result.extend([self.parameters.to_xml_node()])
        return result

    @classmethod
    def from_xml_node(cls, node):
        result = cls(node.get('function'))
        n = get_unique_child(node, ListOfParameters.tag)
        result.parameters = ListOfParameters.from_xml_node(n)
        return result
        
class ListOfEnzymeEfficiencies(ListOf):
    tag = 'listOfEnzymeEfficiencies'
    list_element = EnzymeEfficiency

class TransporterEfficiency(ListOf):
    tag = 'transporterEfficiency'
    list_element = Function
