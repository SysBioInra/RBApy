"""
Module defining macromolecule-specific classes used for RBA XML structures.
"""

# python 2/3 compatiblity
from __future__ import division, print_function, absolute_import

# global imports
from lxml import etree

# local imports
from rba.xml.common import get_unique_child, ListOf

__all__ = ['RbaMacromolecules', 'Component', 'ListOfComponents',
           'Macromolecule', 'ListOfMacromolecules', 'ComponentReference',
           'Composition']


class RbaMacromolecules(object):
    """
    Macromolecule-related structures.

    Attributes
    ----------
    components : ListOfComponents
        List of components of macromolecule (e.g. amino acids for proteins).
    macromolecules : ListOfMacromolecules
        List of macromolecules.
    """

    def __init__(self):
        """
        Default constructor.
        """
        self.components = ListOfComponents()
        self.macromolecules = ListOfMacromolecules()

    def write(self, output_stream, doc_name='RBAMacromolecules'):
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
        root.extend([self.components.to_xml_node(),
                     self.macromolecules.to_xml_node()])
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
        n = get_unique_child(node, ListOfComponents.tag)
        result.components = ListOfComponents.from_xml_node(n)
        n = get_unique_child(node, ListOfMacromolecules.tag)
        result.macromolecules = ListOfMacromolecules.from_xml_node(n)
        return result


class Component(object):
    """
    Component of a macromolecule.

    Attributes
    ----------
    id : str
        Identifier.
    name : str
        Usual name.
    type : str
        Type of molecule (amino acid, vitamin, etc.)
    weight : float
        Weight.
    """

    tag = 'component'

    def __init__(self, id_, name, type_, weight):
        """
        Constructor.

        Parameters
        ----------
        id_ : str
            Identifier.
        name : str
            Usual name.
        type_ : str
            Type of molecule (amino acid, vitamin, etc.)
        weight : float
            Weight.
        """
        self.id = id_
        self.name = name
        self.type = type_
        self.weight = weight

    def to_xml_node(self):
        """
        Convert to xml node.
        """
        result = etree.Element(self.tag)
        result.set('id', self.id)
        result.set('name', self.name)
        result.set('type', self.type)
        result.set('weight', str(self.weight))
        return result

    @classmethod
    def from_xml_node(cls, node):
        """
        Constructor from xml node.
        """
        return cls(node.get('id'), node.get('name'),
                   node.get('type'), node.get('weight'))


class ListOfComponents(ListOf):
    """
    List of Component elements.
    """

    tag = 'listOfComponents'
    list_element = Component


class Macromolecule(object):
    """
    Macromolecule.

    Attributes
    ----------
    id : str
        Identifier.
    compartment : str
        Identifier of compartment where molecule lives.
    composition : Composition
        Composition of macromolecule in terms of components.
    """

    tag = 'macromolecule'

    def __init__(self, id_, compartment, composition=None):
        """
        Constructor.

        Parameters
        ----------
        id_ : str
            Identifier.
        compartment : str
            Identifier of compartment where molecule lives.
        composition : dict, optional
            Dictionary where keys are ids of components and values are their
            stoichiometry within the molecule.
        """
        self.id = id_
        self.compartment = compartment
        self.composition = Composition()
        if composition:
            for comp, sto in composition.items():
                self.composition.append(ComponentReference(comp, sto))

    def to_xml_node(self):
        """
        Convert to xml node.
        """
        result = etree.Element(self.tag)
        result.set('id', self.id)
        result.set('compartment', self.compartment)
        if not self.composition.is_empty():
            result.extend([self.composition.to_xml_node()])
        return result

    @classmethod
    def from_xml_node(cls, node):
        """
        Constructor from xml node.
        """
        result = cls(node.get('id'), node.get('compartment'))
        comp_node = get_unique_child(node, Composition.tag)
        result.composition = Composition.from_xml_node(comp_node)
        return result


class ListOfMacromolecules(ListOf):
    """
    List of Macromolecule elements.
    """

    tag = 'listOfMacromolecules'
    list_element = Macromolecule


class ComponentReference(object):
    """
    Reference to a component, including stoichiometry.

    Attributes
    ----------
    component : str
        Identifier of component.
    stoichiometry : stoichiometry
        Stoichiometry of component.
    """

    tag = 'componentReference'

    def __init__(self, component, stoichiometry):
        """
        Constructor.

        Parameters
        ----------
        component : str
            Identifier of component.
        stoichiometry : stoichiometry
            Stoichiometry of component.
        """
        self.component = component
        self.stoichiometry = stoichiometry

    def to_xml_node(self):
        """
        Convert to xml node.
        """
        result = etree.Element(self.tag)
        result.set('component', self.component)
        result.set('stoichiometry', str(self.stoichiometry))
        return result

    @classmethod
    def from_xml_node(cls, node):
        """
        Constructor from xml node.
        """
        return cls(node.get('component'), float(node.get('stoichiometry')))


class Composition(ListOf):
    """
    List of ComponentReference elements.
    """

    tag = 'composition'
    list_element = ComponentReference
