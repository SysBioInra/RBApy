"""Module defining parameter-specific classes for RBA XML structures."""

# python 2/3 compatiblity
from __future__ import division, print_function, absolute_import

# global imports
from lxml import etree

# local imports
from rba.xml.common import get_unique_child, ListOf

__all__ = ['RbaParameters', 'Parameter', 'ListOfParameters', 'Function',
           'ListOfFunctions', 'FunctionReference', 'ListOfFunctionReferences',
           'Aggregate', 'ListOfAggregates']


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
        """Constructor."""
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
        root.extend([self.functions.to_xml_node(),
                     self.aggregates.to_xml_node()])
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
        n = get_unique_child(node, ListOfFunctions.tag)
        result.functions = ListOfFunctions.from_xml_node(n)
        n = get_unique_child(node, ListOfAggregates.tag)
        result.aggregates = ListOfAggregates.from_xml_node(n)
        return result


class Parameter(object):
    """
    Parameter represented with an id, value couple.

    Attributes
    ----------
    id : str
        Identifier of parameter.
    value : float
        Value of parameter.

    """

    tag = 'parameter'

    def __init__(self, id_, value):
        """
        Constructor.

        Parameters
        ----------
        id_ : str
            Identifier of parameter.
        value : float
            Value of parameter.

        """
        self.id = id_
        self.value = value

    def to_xml_node(self):
        """Convert to xml node."""
        result = etree.Element(self.tag)
        result.set('id', self.id)
        result.set('value', str(self.value))
        return result

    @classmethod
    def from_xml_node(cls, node):
        """Build object from xml node."""
        return cls(node.get('id'), float(node.get('value')))


class ListOfParameters(ListOf):
    """List of Parameter elements."""

    tag = 'listOfParameters'
    list_element = Parameter


class Function(object):
    """
    Function defined by a type and parameters.

    Attributes
    ----------
    id : str or None
        Identifier of function (if applicable).
    type : str
        Type of function.
    parameters : ListOfParameters
        List of parameters used by function.
    variable : str
        Name of variable (if applicable).

    """

    tag = 'function'

    def __init__(self, id_, type_, parameters=None, variable=None):
        """
        Constructor.

        Parameters
        ----------
        id_ : str or None
            identifier of function (if applicable).
        type_: str
            type of function.
        parameters : dict, optional
            dict containing parameters used by function.
        variable : str, optional
            name of variable (set to growth_rate if empty).

        """
        self.id = id_ if id_ is not None else ''
        self.type = type_
        if variable:
            self.variable = variable
        else:
            self.variable = 'growth_rate'
        self.set_parameters(parameters)

    def set_parameters(self, parameters):
        """Create parameter list from dictionary."""
        self.parameters = ListOfParameters()
        if parameters:
            for key, value in parameters.items():
                self.parameters.append(Parameter(key, value))

    def to_xml_node(self):
        """Convert to xml node."""
        result = etree.Element(self.tag)
        result.set('id', self.id)
        result.set('type', self.type)
        if self.variable:
            result.set('variable', self.variable)
        if not self.parameters.is_empty():
            result.extend([self.parameters.to_xml_node()])
        return result

    @classmethod
    def from_xml_node(cls, node):
        """Build object from xml node."""
        result = cls(node.get('id'), node.get('type'),
                     {}, node.get('variable'))
        n = get_unique_child(node, 'listOfParameters', False)
        if n is not None:
            result.parameters = ListOfParameters.from_xml_node(n)
        return result


class ListOfFunctions(ListOf):
    """List of Function elements."""

    tag = 'listOfFunctions'
    list_element = Function


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
        """Convert to xml node."""
        result = etree.Element(self.tag)
        result.set('function', self.function)
        return result

    @classmethod
    def from_xml_node(cls, node):
        """Build object from xml node."""
        return cls(node.get('function'))


class ListOfFunctionReferences(ListOf):
    """List of FunctionReference elements."""

    tag = 'listOfFunctionReferences'
    list_element = FunctionReference


class Aggregate(object):
    """
    Aggregate (composition of Functions).

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
        """Convert to xml node."""
        result = etree.Element(self.tag)
        result.set('id', self.id)
        result.set('type', self.type)
        result.extend([self.function_references.to_xml_node()])
        return result

    @classmethod
    def from_xml_node(cls, node):
        """Build object from xml node."""
        result = cls(node.get('id'), node.get('type'))
        n = get_unique_child(node, ListOfFunctionReferences.tag)
        result.function_references = ListOfFunctionReferences.from_xml_node(n)
        return result


class ListOfAggregates(ListOf):
    """List of Aggregate elements."""

    tag = 'listOfAggregates'
    list_element = Aggregate
