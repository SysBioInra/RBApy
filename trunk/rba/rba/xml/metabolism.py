"""
Module defining metabolism-specific classes used for RBA XML structures.
"""

# python 2/3 compatiblity
from __future__ import division, print_function

# global imports
from lxml import etree

# local imports
from rba.xml.common import (is_true, get_unique_child,
                            ListOf, ListOfProducts, ListOfReactants)


class RbaMetabolism(object):
    """
    Class containing all metabolism-related structures.

    Attributes:
        compartments (ListOfCompartments): list of cell compartments.
        species (ListOfSpecies): list of metabolites.
        reactions (ListOfReactions): list of reactions.
    """
    def __init__(self):
        """
        Default constructor.
        """
        self.compartments = ListOfCompartments()
        self.species = ListOfSpecies()
        self.reactions = ListOfReactions()

    def write(self, output_stream, doc_name = 'RBAMetabolism'):
        """
        Write information as an XML structure.

        Args:
            output_stream: file or buffer where XML structure should be written.
            doc_name: name of XML document.
        """
        root = etree.Element(doc_name)
        root.extend([self.compartments.to_xml_node(),
                     self.species.to_xml_node(), self.reactions.to_xml_node()])
        etree.ElementTree(root).write(output_stream, pretty_print=True)

    @classmethod
    def from_file(cls, input_stream):
        """
        Constructor from XML structure.

        Args:
            input_stream: file or buffer containing XML structure.
        """
        node = etree.ElementTree(file=input_stream).getroot()
        result = cls()
        n = get_unique_child(node, ListOfCompartments.tag)
        result.compartments = ListOfCompartments.from_xml_node(n)
        n = get_unique_child(node, ListOfSpecies.tag)
        result.species = ListOfSpecies.from_xml_node(n)
        n = get_unique_child(node, ListOfReactions.tag)
        result.reactions = ListOfReactions.from_xml_node(n)
        return result


class Compartment(object):
    """
    Class storing compartment information.

    Attributes:
        id: identifier of compartment.
    """

    tag = 'compartment'

    def __init__(self, id_):
        """
        Constructor.

        Args:
            id: identifier of compartment.
        """
        self.id = id_

    def to_xml_node(self):
        """
        Convert to xml node
        """
        result = etree.Element(self.tag)
        result.set('id', self.id)
        return result

    @classmethod
    def from_xml_node(cls, node):
        """
        Constructor from xml node.
        """
        return cls(node.get('id'))
    

class ListOfCompartments(ListOf):
    tag = 'listOfCompartments'
    list_element = Compartment


class Species(object):
    tag = 'species'
    
    def __init__(self, id_, boundary_condition):
        self.id = id_
        self.boundary_condition = boundary_condition

    def to_xml_node(self):
        result = etree.Element(self.tag)
        result.set('id', self.id)
        result.set('boundaryCondition', str(self.boundary_condition).lower())
        return result

    @classmethod
    def from_xml_node(cls, node):
        return cls(node.get('id'), is_true(node.get('boundaryCondition')))


class ListOfSpecies(ListOf):
    tag = 'listOfSpecies'
    list_element = Species


class Reaction(object):
    tag = 'reaction'
    
    def __init__(self, id_, reversible):
        self.id = id_
        self.reversible = reversible
        self.reactants = ListOfReactants()
        self.products = ListOfProducts()

    def to_xml_node(self):
        result = etree.Element(self.tag)
        result.set('id', self.id)
        result.set('reversible', str(self.reversible).lower())
        result.extend([self.reactants.to_xml_node(),
                       self.products.to_xml_node()])
        return result

    @classmethod
    def from_xml_node(cls, node):
        result = cls(node.get('id'), is_true(node.get('reversible')))
        n = get_unique_child(node, ListOfReactants.tag, False)
        if n is not None: result.reactants = ListOfReactants.from_xml_node(n)
        n = get_unique_child(node, ListOfProducts.tag, False)
        if n is not None: result.products = ListOfProducts.from_xml_node(n)
        return result


class ListOfReactions(ListOf):
    tag = 'listOfReactions'
    list_element = Reaction
