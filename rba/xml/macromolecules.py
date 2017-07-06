"""
Module defining macromolecule-specific classes used for RBA XML structures.
"""

# python 2/3 compatiblity
from __future__ import division, print_function

# global imports
from lxml import etree

# local imports
from rba.xml.common import get_unique_child, ListOf

__all__ = ['RbaMacromolecules', 'Component', 'ListOfComponents',
           'Macromolecule', 'ListOfMacromolecules', 'ComponentReference',
           'Composition']


class RbaMacromolecules(object):
    def __init__(self):
        self.components = ListOfComponents()
        self.macromolecules = ListOfMacromolecules()

    def write(self, output_stream, doc_name = 'RBAMacromolecules'):
        root = etree.Element(doc_name)
        root.extend([self.components.to_xml_node(),
                     self.macromolecules.to_xml_node()])
        etree.ElementTree(root).write(output_stream, pretty_print=True)

    @classmethod
    def from_file(cls, input_stream):
        node = etree.ElementTree(file=input_stream).getroot()
        result = cls()
        n = get_unique_child(node, ListOfComponents.tag)
        result.components = ListOfComponents.from_xml_node(n)
        n = get_unique_child(node, ListOfMacromolecules.tag)
        result.macromolecules = ListOfMacromolecules.from_xml_node(n)
        return result

class Component(object):
    tag = 'component'
    
    def __init__(self, id_, name, type_, weight):
        self.id = id_
        self.name = name
        self.type = type_
        self.weight = weight

    def to_xml_node(self):
        result = etree.Element(self.tag)
        result.set('id', self.id)
        result.set('name', self.name)
        result.set('type', self.type)
        result.set('weight', str(self.weight))
        return result

    @classmethod
    def from_xml_node(cls, node):
        return cls(node.get('id'), node.get('name'),
                   node.get('type'), node.get('weight'))

class ListOfComponents(ListOf):
    tag = 'listOfComponents'
    list_element = Component

class Macromolecule(object):
    tag = 'macromolecule'
    
    def __init__(self, id_, compartment, composition=None):
        self.id = id_
        self.compartment = compartment
        self.composition = Composition()
        if composition:
            for comp, sto in composition.items():
                self.composition.append(ComponentReference(comp, sto))

    def to_xml_node(self):
        result = etree.Element(self.tag)
        result.set('id', self.id)
        result.set('compartment', self.compartment)
        if not(self.composition.is_empty()):
            result.extend([self.composition.to_xml_node()])
        return result

    @classmethod
    def from_xml_node(cls, node):
        result = cls(node.get('id'), node.get('compartment'))
        comp_node = get_unique_child(node, Composition.tag)
        result.composition = Composition.from_xml_node(comp_node)
        return result

class ListOfMacromolecules(ListOf):
    tag = 'listOfMacromolecules'
    list_element = Macromolecule

class ComponentReference(object):
    tag = 'componentReference'
    
    def __init__(self, component, stoichiometry):
        self.component = component
        self.stoichiometry = stoichiometry

    def to_xml_node(self):
        result = etree.Element(self.tag)
        result.set('component', self.component)
        result.set('stoichiometry', str(self.stoichiometry))
        return result

    @classmethod
    def from_xml_node(cls, node):
        return cls(node.get('component'), float(node.get('stoichiometry')))

class Composition(ListOf):
    tag = 'composition'
    list_element = ComponentReference
