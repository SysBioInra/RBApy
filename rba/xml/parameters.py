"""
Module defining parameter-specific classes used for RBA XML structures.
"""

# python 2/3 compatiblity
from __future__ import division, print_function

# global imports
from lxml import etree

# local imports
from rba.xml.common import (get_unique_child, ListOf,
                            TargetValue, ListOfFunctions, Parameter)


class RbaParameters(object):
    def __init__(self):
        self.target_densities = ListOfTargetDensities()
        self.functions = ListOfFunctions()
        self.aggregates = ListOfAggregates()

    def write(self, output_stream, doc_name = 'RBAParameters'):
        root = etree.Element(doc_name)
        root.extend([self.target_densities.to_xml_node(),
                     self.functions.to_xml_node(),
                     self.aggregates.to_xml_node()])
        etree.ElementTree(root).write(output_stream, pretty_print=True)

    @classmethod
    def from_file(cls, input_stream):
        node = etree.ElementTree(file=input_stream).getroot()
        result = cls()
        n = get_unique_child(node, ListOfTargetDensities.tag)
        result.target_densities = ListOfTargetDensities.from_xml_node(n)
        n = get_unique_child(node, ListOfFunctions.tag)
        result.functions = ListOfFunctions.from_xml_node(n)
        n = get_unique_child(node, ListOfAggregates.tag)
        result.aggregates = ListOfAggregates.from_xml_node(n)
        return result
    
class TargetDensity(TargetValue):
    tag = 'targetDensity'
    
    def __init__(self, compartment):
        super(TargetDensity, self).__init__()
        self.compartment = compartment

    def to_xml_node(self):
        result = super(TargetDensity,self).to_xml_node()
        result.set('compartment', self.compartment)
        return result
    
    @classmethod
    def from_xml_node(cls, node):
        result = cls(node.get('compartment'))
        result._init_from_xml_node(node)
        return result

class ListOfTargetDensities(ListOf):
    tag = 'listOfTargetDensities'
    list_element = TargetDensity

class FunctionReference(object):
    tag = 'functionReference'

    def __init__(self, function):
        self.function = function

    def to_xml_node(self):
        result = etree.Element(self.tag)
        result.set('function', self.function)
        return result

    @classmethod
    def from_xml_node(cls, node):
        return cls(node.get('function'))

class ListOfFunctionReferences(ListOf):
    tag = 'listOfFunctionReferences'
    list_element = FunctionReference
    
class Aggregate(object):
    tag = 'aggregate'

    def __init__(self, id_, type_):
        if id_ is not None: self.id = id_
        else: self.id = ''
        self.type = type_
        self.function_references = ListOfFunctionReferences()

    def to_xml_node(self):
        result = etree.Element(self.tag)
        result.set('id', self.id)
        result.set('type', self.type)
        result.extend([self.function_references.to_xml_node()])
        return result

    @classmethod
    def from_xml_node(cls, node):
        result = cls(node.get('id'), node.get('type'))
        n = get_unique_child(node, ListOfFunctionReferences.tag)
        result.function_references = ListOfFunctionReferences.from_xml_node(n)
        return result

class ListOfAggregates(ListOf):
    tag = 'listOfAggregates'
    list_element = Aggregate
