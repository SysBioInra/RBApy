"""
Module defining parameter-specific classes used for RBA XML structures.
"""

# python 2/3 compatiblity
from __future__ import division, print_function

# global imports
from lxml import etree

# local imports
from rba.xml.common import (get_unique_child, ListOf,
                            TargetValue, ListOfFunctions)

__all__ = ['RbaParameters', 'FunctionReference', 'ListOfFunctionReferences',
           'TargetDensity', 'ListOfTargetDensities', 'Aggregate',
           'ListOfAggregates']


class RbaParameters(object):
    """
    Parameter-related structures.
    
    Attributes
    ----------
    target_densities : ListOfTargetDensities
        List of density constraints.
    functions : ListOfFunctions
        List of user-defined functions.
    aggregates : ListOfAggregated
        List of user-defined aggregates (composition of functions).
    """

    def __init__(self):
        """
        Default constructor.
        """
        self.target_densities = ListOfTargetDensities()
        self.functions = ListOfFunctions()
        self.aggregates = ListOfAggregates()

    def write(self, output_stream, doc_name='RBAParameters'):
        """
        Write information as an XML structure.

        Parameters
        ----------
        output_stream : file or buffer
            Location where XML structure should be written.
        doc_name : str, optional
            Name of XML document.
        """
        root = etree.Element(doc_name)
        root.extend([self.target_densities.to_xml_node(),
                     self.functions.to_xml_node(),
                     self.aggregates.to_xml_node()])
        etree.ElementTree(root).write(output_stream, pretty_print=True)

    @classmethod
    def from_file(cls, input_stream):
        """
        Constructor from XML structure.

        Parameters
        ----------
        input_stream : file or buffer
            Location containing XML structure.
        """
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
    """
    Density constraint.

    Attributes
    ----------
    compartment : str
        Compartment identifier.
    """

    tag = 'targetDensity'

    def __init__(self, compartment):
        """
        Constructor.

        Parameters
        ----------
        compartment : str
            Compartment identifier.
        """
        super(TargetDensity, self).__init__()
        self.compartment = compartment

    def to_xml_node(self):
        """
        Convert to xml node.
        """
        result = super(TargetDensity, self).to_xml_node()
        result.set('compartment', self.compartment)
        return result

    @classmethod
    def from_xml_node(cls, node):
        """
        Constructor from xml node.
        """
        result = cls(node.get('compartment'))
        result._init_from_xml_node(node)
        return result


class ListOfTargetDensities(ListOf):
    """
    List of TargetDensity elements.
    """

    tag = 'listOfTargetDensities'
    list_element = TargetDensity


class FunctionReference(object):
    """
    Reference to a Function.

    Attributes
    ----------
    function : str
        Function identifier.
    """

    tag = 'functionReference'

    def __init__(self, function):
        """
        Constructor.

        Parameters
        ----------
        function : str
            Function identifier.
        """
        self.function = function

    def to_xml_node(self):
        """
        Convert to xml node.
        """
        result = etree.Element(self.tag)
        result.set('function', self.function)
        return result

    @classmethod
    def from_xml_node(cls, node):
        """
        Constructor from xml node.
        """
        return cls(node.get('function'))


class ListOfFunctionReferences(ListOf):
    """
    List of FunctionReference elements.
    """

    tag = 'listOfFunctionReferences'
    list_element = FunctionReference


class Aggregate(object):
    """
    Aggregate (composition of Functions)

    Attributes
    ----------
    id : str
        Identifier.
    type : str
        Type of aggregation (e.g. 'multiplication').
    function_references : ListOfFunctionReferences
        References to functions aggregated.
    """

    tag = 'aggregate'

    def __init__(self, id_, type_):
        """
        Constructor.

        Parameters
        ----------
        id_ : str
            Identifier.
        type_ : str
            Type of aggregation (e.g. 'multiplication').
        """
        self.id = id_ if id_ is not None else ''
        self.type = type_
        self.function_references = ListOfFunctionReferences()

    def to_xml_node(self):
        """
        Convert to xml node.
        """
        result = etree.Element(self.tag)
        result.set('id', self.id)
        result.set('type', self.type)
        result.extend([self.function_references.to_xml_node()])
        return result

    @classmethod
    def from_xml_node(cls, node):
        """
        Constructor from xml node.
        """
        result = cls(node.get('id'), node.get('type'))
        n = get_unique_child(node, ListOfFunctionReferences.tag)
        result.function_references = ListOfFunctionReferences.from_xml_node(n)
        return result


class ListOfAggregates(ListOf):
    """
    List of Aggregate elements.
    """

    tag = 'listOfAggregates'
    list_element = Aggregate
