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
    """
    Enzyme-related structures.

    Attributes
    ----------
    efficiency_functions : ListOfEfficiencyFunctions
        Efficiency function families for enzymes.
    enzymes : ListOfEnzymes
        Enzymes in the system.
    """

    def __init__(self):
        """
        Default constructor.
        """
        self.efficiency_functions = ListOfEfficiencyFunctions()
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
        root.extend([self.efficiency_functions.to_xml_node(),
                     self.enzymes.to_xml_node()])
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
        n = get_unique_child(node, ListOfEfficiencyFunctions.tag)
        result.efficiency_functions = ListOfEfficiencyFunctions.from_xml_node(n)
        n = get_unique_child(node, ListOfEnzymes.tag)
        result.enzymes = ListOfEnzymes.from_xml_node(n)
        return result


class ListOfEfficiencyFunctions(ListOf):
    """
    List of Function elements representing efficiency function families.
    """

    tag = 'listOfEfficiencyFunctions'
    list_element = Function


class Enzyme(object):
    """
    Enzyme.

    Attributes
    ----------
    zero_cost : bool
        True if enzyme should be synthesized for free.
    id : str
        Identifier.
    machinery_composition : MachineryComposition
        Composition of enzyme in terms of proteins.
    enzymatic_activity : EnzymaticActivity
        Enzymatic parameters (one set of parameters for each efficiency
        function family defined) and transporter activity (if applicable).
    """

    tag = 'enzyme'

    def __init__(self, id_, zero_cost=False):
        """
        Constructor.

        Parameters
        ----------
        id_ : str
            Identifier.
        zero_cost : bool, optional
            True if enzyme should be synthesized for free.
        """
        self.id = id_
        self.zero_cost = zero_cost if zero_cost is not None else False
        self.machinery_composition = MachineryComposition()
        self.enzymatic_activity = EnzymaticActivity(id_)

    def to_xml_node(self):
        """
        Convert to xml node.
        """
        result = etree.Element(self.tag)
        result.set('id', self.id)
        result.set('zeroCost', str(self.zero_cost).lower())
        if not self.machinery_composition.is_empty():
            result.extend([self.machinery_composition.to_xml_node()])
        result.extend([self.enzymatic_activity.to_xml_node()])
        return result

    @classmethod
    def from_xml_node(cls, node):
        """
        Constructor from xml node.
        """
        result = cls(node.get('id'), is_true(node.get('zeroCost')))
        n = get_unique_child(node, MachineryComposition.tag, False)
        if n is not None:
            result.machinery_composition = MachineryComposition.from_xml_node(n)
        n = get_unique_child(node, EnzymaticActivity.tag)
        result.enzymatic_activity = EnzymaticActivity.from_xml_node(n)
        return result


class ListOfEnzymes(ListOf):
    """
    List of Enzyme elements.
    """

    tag = 'listOfEnzymes'
    list_element = Enzyme


class EnzymaticActivity(object):
    """
    Enzymatic parameters for efficiency families and transporting.

    Attributes
    ----------
    reaction : str
        Identifier of reaction catalyzed.
    enzyme_efficiencies : ListOfEnzymeEfficiencies
        List of enzymatic parameters (one set for every efficiency function).
    transporter_efficiency : TransporterEfficiency
        Parameter for transporting activity (if applicable).
    """

    tag = 'enzymaticActivity'

    def __init__(self, reaction):
        """
        Constructor.

        Parameters
        ----------
        reaction : str
            Identifier of reaction catalyzed.
        """
        self.reaction = reaction
        self.enzyme_efficiencies = ListOfEnzymeEfficiencies()
        self.transporter_efficiency = TransporterEfficiency()

    def to_xml_node(self):
        """
        Convert to xml node.
        """
        result = etree.Element(self.tag)
        result.set('reaction', self.reaction)
        # optional subelements
        opt = [self.enzyme_efficiencies, self.transporter_efficiency]
        result.extend([e.to_xml_node() for e in opt if not e.is_empty()])
        return result

    @classmethod
    def from_xml_node(cls, node):
        """
        Constructor from xml node.
        """
        reaction = node.get('reaction')
        result = cls(reaction)
        n = get_unique_child(node, ListOfEnzymeEfficiencies.tag)
        result.enzyme_efficiencies = ListOfEnzymeEfficiencies.from_xml_node(n)
        n = get_unique_child(node, TransporterEfficiency.tag, False)
        if n is not None:
            result.transporter_efficiency = (TransporterEfficiency
                                             .from_xml_node(n))
        return result


class EnzymeEfficiency(object):
    """
    Enzymatic parameter of a given efficiency function.

    Attributes
    ----------
    function : str
        Identifier of efficiency function for which parameters apply.
    parameters : ListOfParameters
        Parameters used to compute efficiency function for given enzyme.
    """

    tag = 'enzymeEfficiency'

    def __init__(self, function, parameters=None):
        """
        Constructor.

        Parameters
        ----------
        function : str
            Identifier of efficiency function for which parameters apply.
        parameters : dict, optional
            Dictionary mapping parameter ids with their value.
        """
        self.function = function
        self.parameters = ListOfParameters()
        if parameters:
            for key, value in parameters.items():
                self.parameters.append(Parameter(key, value))

    def to_xml_node(self):
        """
        Convert to xml node.
        """
        result = etree.Element(self.tag)
        result.set('function', self.function)
        result.extend([self.parameters.to_xml_node()])
        return result

    @classmethod
    def from_xml_node(cls, node):
        """
        Constructor from xml node.
        """
        result = cls(node.get('function'))
        n = get_unique_child(node, ListOfParameters.tag)
        result.parameters = ListOfParameters.from_xml_node(n)
        return result


class ListOfEnzymeEfficiencies(ListOf):
    """
    List of EnzymeEfficiency elements.
    """

    tag = 'listOfEnzymeEfficiencies'
    list_element = EnzymeEfficiency


class TransporterEfficiency(ListOf):
    """
    List of Function elements used to define transporter efficiency.
    """

    tag = 'transporterEfficiency'
    list_element = Function
