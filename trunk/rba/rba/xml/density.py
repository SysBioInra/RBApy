"""Module defining density-specific classes for RBA XML structures."""

# python 2/3 compatiblity
from __future__ import division, print_function, absolute_import

# global imports
from lxml import etree

# local imports
from rba.xml.common import get_unique_child, ListOf, TargetValue

__all__ = ['RbaDensity', 'TargetDensity', 'ListOfTargetDensities']


class RbaDensity(object):
    """
    Parameter-related structures.

    Attributes
    ----------
    target_densities : ListOfTargetDensities
        List of density constraints.

    """

    def __init__(self):
        """Constructor."""
        self.target_densities = ListOfTargetDensities()

    def write(self, output_stream, doc_name='RBADensity'):
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
        root.extend([self.target_densities.to_xml_node()])
        etree.ElementTree(root).write(output_stream, pretty_print=True)

    @classmethod
    def from_file(cls, input_stream):
        """
        Build object from XML structure.

        Parameters
        ----------
        input_stream : file or buffer
            Location containing XML structure.

        """
        node = etree.ElementTree(file=input_stream).getroot()
        result = cls()
        n = get_unique_child(node, ListOfTargetDensities.tag)
        result.target_densities = ListOfTargetDensities.from_xml_node(n)
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
        """Convert to xml node."""
        result = super(TargetDensity, self).to_xml_node()
        result.set('compartment', self.compartment)
        return result

    @classmethod
    def from_xml_node(cls, node):
        """Build object from xml node."""
        result = cls(node.get('compartment'))
        result._init_from_xml_node(node)
        return result


class ListOfTargetDensities(ListOf):
    """List of TargetDensity elements."""

    tag = 'listOfTargetDensities'
    list_element = TargetDensity
