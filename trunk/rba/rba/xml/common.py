"""
Module defining common function and classes used for RBA XML structures.
"""

# python 2/3 compatiblity
from __future__ import division, print_function

# global imports
from lxml import etree


def is_true(attribute):
    """
    Determine whether an XML attribute evaluates as the logical true value.

    Args:
        attribute: attribute value.

    Returns:
        True if attribute value is considered to represent the 'true' boolean
        state, False otherwise.
    """
    return attribute.lower() == 'true' or attribute == '1'


def get_unique_child(parent, child_name, strict=True):
    """
    Get children node matching given name.

    Args:
        parent: Parent node.
        child_name: Name of the node to look for.
        strict: Boolean flag telling whether there should be strictly one
            child. If set to True, an error is raised if no child was found.

    Returns:
        Child if exactly one node is found. If no child is found,
        None is returned if strict was set to False.

    Raises:
        If no child is found, a UserWarning exception is raised if strict 
        is set to True.
        A UserWarning exception is raised if more than one node is found.
    """
    children = parent.findall(child_name)
    if len(children) == 1:
        return children[0]
    elif not children and not strict:
        return None
    else:
        if strict:
            raise UserWarning('A {} should have exactly one {} node child.'
                              .format(parent.tag, child_name))
        else:
            raise UserWarning('A {} should not have more than one {} node '
                              'children.'.format(parent.tag, child_name))


class ListOf(object):
    """
    Abstract class used to store list of XML objects of same type.
    
    Attributes:
        list_element: type of elements stored by list.
    """
    
    list_element = None
    
    def __init__(self):
        """
        Default constructor.
        """
        self._elements = []

    def __getitem__(self, i):
        """
        Return item with given index.
        """
        return self._elements[i]

    def __len__(self):
        """
        Return length of list.
        """
        return len(self._elements)

    def append(self, element):
        """
        Append element to list.
        """
        assert(isinstance(element, self.list_element))
        self._elements.append(element)

    def is_empty(self):
        """
        Return whether list is empty.
        """
        return len(self._elements) == 0

    def to_xml_node(self):
        """
        Convert to xml node.
        """
        result = etree.Element(self.tag)
        result.extend([e.to_xml_node() for e in self._elements])
        return result

    @classmethod
    def from_xml_node(cls, node):
        """
        Constructor from xml node.
        """
        result = cls()
        for n in node.iterfind(cls.list_element.tag):
            result.append(cls.list_element.from_xml_node(n))
        return result


class MachineryComposition(object):
    """
    Class storing machinery composition information.

    Attributes:
        reactants (ListOfReactants): list of reactants used to assemble
            machinery.
        products (ListOfProducts): list of byproducts generated while
            assembling machinery.
    """
    
    tag = 'machineryComposition'
    
    def __init__(self):
        """
        Default constructor.
        """
        self.reactants = ListOfReactants()
        self.products = ListOfProducts()

    def is_empty(self):
        """
        Return whether composition is fully empty.
        """
        return self.reactants.is_empty() and self.products.is_empty()

    def to_xml_node(self):
        """
        Convert to xml node.
        """
        result = etree.Element(self.tag)
        # optional subelements
        opt = [self.reactants, self.products]
        result.extend([e.to_xml_node() for e in opt if not(e.is_empty())])
        return result

    @classmethod
    def from_xml_node(cls, node):
        """
        Constructor from xml node.
        """
        result = cls()
        n = get_unique_child(node, ListOfReactants.tag, False)
        if n is not None:
            result.reactants = ListOfReactants.from_xml_node(n)
        n = get_unique_child(node, ListOfProducts.tag, False)
        if n is not None:
            result.products = ListOfProducts.from_xml_node(n)
        return result
    

class SpeciesReference(object):
    """
    Class storing reference to a species and its stoichiometry.

    Attributes:
        species (str): identifier of species.
        stoichiometry (float): stoichiometry of species.
    """
    
    tag = 'speciesReference'
    
    def __init__(self, species, stoichiometry):
        """
        Constructor.

        Args:
            species (str): identifier of species.
            stoichiometry (float): stoichiometry of species.        
        """
        self.species = species
        self.stoichiometry = stoichiometry

    def to_xml_node(self):
        """
        Convert to xml node.
        """
        result = etree.Element(self.tag)
        result.set('species', self.species)
        result.set('stoichiometry', str(self.stoichiometry))
        return result

    @classmethod
    def from_xml_node(cls, node):
        """
        Constructor from xml node.
        """
        return cls(node.get('species'), float(node.get('stoichiometry')))


class ListOfReactants(ListOf):
    """
    Class storing a list of SpeciesReference representing reactants.
    """
    
    tag = 'listOfReactants'
    list_element = SpeciesReference


class ListOfProducts(ListOf):
    """
    Class storing a list of SpeciesReference representing products.
    """
    
    tag = 'listOfProducts'
    list_element = SpeciesReference


class Parameter(object):
    """
    Class storing parameter values.

    Attributes:
        id: identifier of parameter.
        value: value of parameter.
    """
    
    tag = 'parameter'
    
    def __init__(self, id_, value):
        """
        Constructor.
        
        Args:
            id_: identifier of parameter.
            value: value of parameter.
        """
        self.id = id_
        self.value = value

    def to_xml_node(self):
        """
        Convert to xml node.
        """
        result = etree.Element(self.tag)
        result.set('id', self.id)
        result.set('value', str(self.value))
        return result

    @classmethod
    def from_xml_node(cls, node):
        """
        Constructor from xml node.
        """
        return cls(node.get('id'), float(node.get('value')))


class ListOfParameters(ListOf):
    """
    Class storing a list of Parameters.
    """
    
    tag = 'listOfParameters'
    list_element = Parameter


class Function(object):
    """
    Class storing a Function.

    Attributes:
        id: identifier of function (if applicable).
        type: type of function.
        parameters (ListOfParameters): list of parameters used by function.
        variable: name of variable (if applicable).
    """
    
    tag = 'function'
    
    def __init__(self, id_, type_, parameters=None, variable=''):
        """
        Constructor.
        
        Attributes:
            id: identifier of function (if applicable).
            type: type of function.
            parameters (ListOfParameters): list of parameters used by function.
            variable: name of variable (if applicable).
        """
        if id_ is not None: self.id = id_
        else: self.id = ''
        self.type = type_
        self.variable = variable
        self.parameters = ListOfParameters()
        if parameters:
            for key, value in parameters.items():
                self.parameters.append(Parameter(key, value))

    def to_xml_node(self):
        """
        Convert to xml node.
        """
        result = etree.Element(self.tag)
        result.set('id', self.id)
        result.set('type', self.type)
        if self.variable: result.set('variable', self.variable)
        if not(self.parameters.is_empty()):
            result.extend([self.parameters.to_xml_node()])
        return result

    @classmethod
    def from_xml_node(cls, node):
        """
        Constructor from xml node.
        """
        result = cls(node.get('id'), node.get('type'), {}, node.get('variable'))
        n = get_unique_child(node, 'listOfParameters', False)
        if n is not None:
            result.parameters = ListOfParameters.from_xml_node(n)
        return result


class ListOfFunctions(ListOf):
    """
    Class storing a list of Functions.
    """
    
    tag = 'listOfFunctions'
    list_element = Function


class TargetValue(object):
    """
    Class storing a target value.

    Attributes:
        value: exact target value (None if no exact value to match).
        lower_bound: lower bound on target value (None if no lower bound).
        upper_bound: upper bound on target value (None if no upper bound).
    """
    tag = 'targetValue'
    
    def __init__(self):
        """
        Default constructor.
        """
        self.value = None
        self.lower_bound = None
        self.upper_bound = None

    def is_empty(self):
        """
        Return whether all attributes are unspecified.
        """
        return self.value is None and self.lower_bound is None \
            and self.upper_bound is None

    def to_xml_node(self):
        """
        Convert to xml node.
        """
        result = etree.Element(self.tag)
        if self.value is not None:
            result.set('value', str(self.value))
        if self.lower_bound is not None:
            result.set('lowerBound', str(self.lower_bound))
        if self.upper_bound is not None:
            result.set('upperBound', str(self.upper_bound))
        return result
    
    @classmethod
    def from_xml_node(cls, node):
        """
        Constructor from xml node.
        """
        result = cls()
        result._init_from_xml_node(node)
        return result

    def _init_from_xml_node(self, node):
        """
        Match attributes with given node.
        """
        self.value = node.get('value')
        self.lower_bound = node.get('lowerBound')
        self.upper_bound = node.get('upperBound')
