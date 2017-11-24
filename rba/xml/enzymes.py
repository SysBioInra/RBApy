"""Module defining enzyme-specific classes used for RBA XML structures."""

# python 2/3 compatiblity
from __future__ import division, print_function, absolute_import

# global imports
from lxml import etree

# local imports
from rba.xml.common import (
    get_unique_child, is_true, ListOf, MachineryComposition
    )

__all__ = ['RbaEnzymes', 'Enzyme', 'ListOfEnzymes']


class RbaEnzymes(object):
    """
    Enzyme-related structures.

    Attributes
    ----------
    enzymes : ListOfEnzymes
        Enzymes in the system.

    """

    def __init__(self):
        """Build default object."""
        self.enzymes = ListOfEnzymes()

    def write(self, output_stream, doc_name='RBAEnzymes'):
        """
        Write information as an XML structure.

        Parameters
        ----------
        output_stream : file or buffer
            Location where XML structure should be written.
        doc_name : str
            Name of XML document.

        """
        root = etree.Element(doc_name)
        root.extend([self.enzymes.to_xml_node()])
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
        n = get_unique_child(node, ListOfEnzymes.tag)
        result.enzymes = ListOfEnzymes.from_xml_node(n)
        return result


class Enzyme(object):
    """
    Enzyme.

    Attributes
    ----------
    id : str
        Identifier.
    reaction : str
        Identifier of reaction catalyzed.
    forward_efficiency : str
        Identifier of function describing forward efficiency.
    backward_efficiency : str
        Identifier of function describing backward efficiency.
    machinery_composition : MachineryComposition
        Composition of enzyme in terms of proteins.
    zero_cost : bool
        True if enzyme should be synthesized for free.

    """

    tag = 'enzyme'

    def __init__(self, id_, reaction, forward_efficiency,
                 backward_efficiency, zero_cost=False):
        """
        Constructor.

        Parameters
        ----------
        id_ : str
            Identifier.
        reaction : str
            Identifier of reaction catalyzed.
        froward_efficiency : str
            Identifier of function describing forward efficiency.
        backward_efficiency : str
            Identifier of function describing backward efficiency.
        zero_cost : bool
            True if enzyme should be synthesized for free (False by default).

        """
        self.id = id_
        self.reaction = reaction
        self.forward_efficiency = forward_efficiency
        self.backward_efficiency = backward_efficiency
        self.zero_cost = zero_cost if zero_cost is not None else False
        self.machinery_composition = MachineryComposition()

    def to_xml_node(self):
        """Convert to xml node."""
        result = etree.Element(self.tag)
        result.set('id', self.id)
        result.set('reaction', self.reaction)
        result.set('forward_efficiency', self.forward_efficiency)
        result.set('backward_efficiency', self.backward_efficiency)
        result.set('zeroCost', str(self.zero_cost).lower())
        if not self.machinery_composition.is_empty():
            result.extend([self.machinery_composition.to_xml_node()])
        return result

    @classmethod
    def from_xml_node(cls, node):
        """Build from xml node."""
        result = cls(node.get('id'), node.get('reaction'),
                     node.get('forward_efficiency'),
                     node.get('backward_efficiency'),
                     is_true(node.get('zeroCost')))
        n = get_unique_child(node, MachineryComposition.tag, False)
        if n is not None:
            result.machinery_composition \
                = MachineryComposition.from_xml_node(n)
        return result


class ListOfEnzymes(ListOf):
    """List of Enzyme elements."""

    tag = 'listOfEnzymes'
    list_element = Enzyme
