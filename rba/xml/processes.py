"""Module defining process-specific classes used for RBA XML structures."""

# python 2/3 compatiblity
from __future__ import division, print_function, absolute_import

# global imports
from lxml import etree

# local imports
from rba.xml.common import (get_unique_child, ListOf, SpeciesReference,
                            TargetValue, MachineryComposition,
                            ListOfReactants, ListOfProducts)

__all__ = ['RbaProcesses', 'Process', 'ListOfProcesses',
           'Machinery', 'Capacity', 'Processing', 'Processings',
           'ListOfInputs', 'ListOfProductions', 'ListOfDegradations',
           'ProcessingMap', 'ListOfProcessingMaps',
           'ConstantProcessing', 'ComponentProcessing',
           'ListOfComponentProcessings']


class RbaProcesses(object):
    """
    Process-related structures.

    Attributes
    ----------
    processes : ListOfProcesses
        List of cell processes.
    processing_maps : ListOfProcessingMaps
        List of chemical reactions used for macromolecule synthesis.

    """

    def __init__(self):
        """Constructor."""
        self.processes = ListOfProcesses()
        self.processing_maps = ListOfProcessingMaps()

    def write(self, output_stream, doc_name='RBAProcesses'):
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
        root.extend([self.processes.to_xml_node(),
                     self.processing_maps.to_xml_node()])
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
        n = get_unique_child(node, ListOfProcesses.tag)
        result.processes = ListOfProcesses.from_xml_node(n)
        n = get_unique_child(node, ListOfProcessingMaps.tag)
        result.processing_maps = ListOfProcessingMaps.from_xml_node(n)
        return result


class Process(object):
    """
    Cell process.

    Attributes
    ----------
    id : str
        Identifier.
    name : str
        Usual name.
    machinery : Machinery
        Molecular machinery used by process (can be empty).
    processings : rba.xml.Processings
        Molecules produced / degraded by process (can be empty).

    """

    tag = 'process'

    def __init__(self, id_, name):
        """
        Constructor.

        Parameters
        ----------
        id_ : str
            Identifier.
        name : name
            Usual name.

        """
        self.id = id_
        self.name = name
        self.machinery = Machinery()
        self.processings = Processings()

    def to_xml_node(self):
        """Convert to xml node."""
        process = etree.Element(self.tag)
        process.set('id', self.id)
        process.set('name', self.name)
        # optional subelements
        opt = [self.machinery, self.processings]
        process.extend([e.to_xml_node() for e in opt if not e.is_empty()])
        return process

    @classmethod
    def from_xml_node(cls, node):
        """Build object from xml node."""
        result = cls(node.get('id'), node.get('name'))
        n = get_unique_child(node, Machinery.tag, False)
        if n is not None:
            result.machinery = Machinery.from_xml_node(n)
        n = get_unique_child(node, Processings.tag, False)
        if n is not None:
            result.processings = Processings.from_xml_node(n)
        return result


class ListOfProcesses(ListOf):
    """List of Process elements."""

    tag = 'listOfProcesses'
    list_element = Process


class Machinery(object):
    """
    Molecular machinery.

    Attributes
    ----------
    machinery_composition : MachineryComposition
        Molecular composition.
    capacity : Capacity
        Processing capacity (number of units processed per unit of time).

    """

    tag = 'machinery'

    def __init__(self):
        """Constructor."""
        self.machinery_composition = MachineryComposition()
        self.capacity = Capacity()

    def is_empty(self):
        """Return whether moleculer composition is empty."""
        return self.machinery_composition.is_empty() \
            and self.capacity.is_empty()

    def to_xml_node(self):
        """Convert to xml node."""
        result = etree.Element(self.tag)
        result.extend([self.machinery_composition.to_xml_node(),
                       self.capacity.to_xml_node()])
        return result

    @classmethod
    def from_xml_node(cls, node):
        """Build object from xml node."""
        result = cls()
        machinery_node = get_unique_child(node, MachineryComposition.tag)
        result.machinery_composition \
            = MachineryComposition.from_xml_node(machinery_node)
        capacity_node = get_unique_child(node, Capacity.tag)
        result.capacity = Capacity.from_xml_node(capacity_node)
        return result


class Capacity(TargetValue):
    """
    Machinery capacity.

    Capacity is used to write constraints on process machineries.
    If a value is specified, a strict constraint is implemented,
    meaning that all machinery units must work at full capacity.
    If an upper bound is specified, the constraint becomes an inequality.

    """

    tag = 'capacity'


class Processings(object):
    """
    Set of production and degradation operations.

    Attributes
    ----------
    productions : ListOfProductions
        Production operations.
    degradations : ListOfDegradations
        Degradation operations.

    """

    tag = 'processings'

    def __init__(self):
        """Constructor."""
        self.productions = ListOfProductions()
        self.degradations = ListOfDegradations()

    def is_empty(self):
        """Return whether there are no operations."""
        return self.productions.is_empty() and self.degradations.is_empty()

    def to_xml_node(self):
        """Convert to xml node."""
        result = etree.Element(self.tag)
        # optional subelements
        opt = [self.productions, self.degradations]
        result.extend([e.to_xml_node() for e in opt if not e.is_empty()])
        return result

    @classmethod
    def from_xml_node(cls, node):
        """Build object from xml node."""
        result = cls()
        n = get_unique_child(node, ListOfProductions.tag, False)
        if n is not None:
            result.productions = ListOfProductions.from_xml_node(n)
        n = get_unique_child(node, ListOfDegradations.tag, False)
        if n is not None:
            result.degradations = ListOfDegradations.from_xml_node(n)
        return result


class ListOfInputs(ListOf):
    """List of species inputs."""

    tag = 'listOfInputs'
    list_element = SpeciesReference


class Processing(object):
    """
    Reference to a production or degradation operation.

    Attributes
    ----------
    processing_map : str
        Identifier of ProcessingMap defining metabolites involved in operation.
    set : str
        Identifier of Macromolecule set being processed.
    inputs : rba.xml.ListOfInputs
        List of Macromolecules being operated on.

    """

    tag = 'processing'

    def __init__(self, map_, set_):
        """
        Constructor.

        Parameters
        ----------
        map_ :  str
            Identifier of ProcessingMap defining metabolites involved in
            operation.
        set : str
            Identifier of Macromolecule set being processed.

        """
        self.processing_map = map_
        self.set = set_
        self.inputs = ListOfInputs()

    def to_xml_node(self):
        """Convert to xml node."""
        result = etree.Element(self.tag)
        result.set('processingMap', self.processing_map)
        result.set('set', self.set)
        result.extend([self.inputs.to_xml_node()])
        return result

    @classmethod
    def from_xml_node(cls, node):
        """Build object from xml node."""
        result = cls(node.get('processingMap'), node.get('set'))
        n = get_unique_child(node, ListOfInputs.tag)
        result.inputs = ListOfInputs.from_xml_node(n)
        return result


class ListOfProductions(ListOf):
    """List of Operation elements representing productions."""

    tag = 'listOfProductions'
    list_element = Processing


class ListOfDegradations(ListOf):
    """List of Operation elements representing degradations."""

    tag = 'listOfDegradations'
    list_element = Processing


class ProcessingMap(object):
    """
    Cost of [un]polymerization of components of a macromolecule.

    Attributes
    ----------
    id : str
        Identifier of map.
    constant_processing : ConstantProcessing
        Cost associated with synthesis/degradation that do not depend
        on the number of components (e.g. initialization).
    component_processings : ListOfComponentProcessings
        Cost associated with [un]polymerization of components.

    """

    tag = 'processingMap'

    def __init__(self, id_):
        """
        Constructor.

        Parameters
        ----------
        id_ : str
            Identifier of map.

        """
        self.id = id_
        self.constant_processing = ConstantProcessing()
        self.component_processings = ListOfComponentProcessings()

    def to_xml_node(self):
        """Convert to xml node."""
        result = etree.Element(self.tag)
        result.set('id', self.id)
        result.extend([self.constant_processing.to_xml_node(),
                       self.component_processings.to_xml_node()])
        return result

    @classmethod
    def from_xml_node(cls, node):
        """Build object from xml node."""
        result = cls(node.get('id'))
        n = get_unique_child(node, ConstantProcessing.tag, False)
        if n is not None:
            result.constant_processing = ConstantProcessing.from_xml_node(n)
        n = get_unique_child(node, ListOfComponentProcessings.tag)
        result.component_processings \
            = ListOfComponentProcessings.from_xml_node(n)
        return result


class ListOfProcessingMaps(ListOf):
    """List of ProcessingMap elements."""

    tag = 'listOfProcessingMaps'
    list_element = ProcessingMap


class ConstantProcessing(object):
    """
    Metabolites consumed/produced independent of component number/nature.

    Attributes
    ----------
    reactants : ListOfReactants
        Metabolites consumed.
    products : ListOfProducts
        Metabolites generated as byproducts.

    """

    tag = 'constantProcessing'

    def __init__(self):
        """Constructor."""
        self.reactants = ListOfReactants()
        self.products = ListOfProducts()

    def to_xml_node(self):
        """Convert to xml node."""
        result = etree.Element(self.tag)
        # optional subelements
        opt = [self.reactants, self.products]
        result.extend([e.to_xml_node() for e in opt if not e.is_empty()])
        return result

    @classmethod
    def from_xml_node(cls, node):
        """Build object from xml node."""
        result = cls()
        n = get_unique_child(node, ListOfReactants.tag, False)
        if n is not None:
            result.reactants = ListOfReactants.from_xml_node(n)
        n = get_unique_child(node, ListOfProducts.tag, False)
        if n is not None:
            result.products = ListOfProducts.from_xml_node(n)
        return result


class ComponentProcessing(object):
    """
    Cost of [un]polymerization of given component.

    Attributes
    ----------
    component : str
        Identifier of component.
    machinery_cost : float
        Capacity cost for process machinery to [un]polymerize component.
    reactants : ListOfReactants
        Metabolites consumed.
    products : ListOfProducts
        Metabolites generated as byproducts.

    """

    tag = 'componentProcessing'

    def __init__(self, component, machinery_cost=0.0):
        """
        Constructor.

        Parameters
        ----------
        component : str
            Identifier of component.
        machinery_cost : float, optional
            Capacity cost for process machinery to [un]polymerize component.

        """
        self.component = component
        self.machinery_cost = machinery_cost
        self.reactants = ListOfReactants()
        self.products = ListOfProducts()

    def to_xml_node(self):
        """Convert to xml node."""
        result = etree.Element(self.tag)
        result.set('component', self.component)
        result.set('machineryCost', str(self.machinery_cost))
        # optional subelements
        opt = [self.reactants, self.products]
        result.extend([e.to_xml_node() for e in opt if not e.is_empty()])
        return result

    @classmethod
    def from_xml_node(cls, node):
        """Build object from xml node."""
        try:
            machinery_cost = float(node.get('machineryCost'))
        except TypeError:
            machinery_cost = 0.0
        result = cls(node.get('component'), machinery_cost)
        n = get_unique_child(node, ListOfReactants.tag, False)
        if n is not None:
            result.reactants = ListOfReactants.from_xml_node(n)
        n = get_unique_child(node, ListOfProducts.tag, False)
        if n is not None:
            result.products = ListOfProducts.from_xml_node(n)
        return result


class ListOfComponentProcessings(ListOf):
    """List of ComponentProcessing elements."""

    tag = 'listOfComponentProcessings'
    list_element = ComponentProcessing
