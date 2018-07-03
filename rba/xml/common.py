"""Common function and classes used for RBA XML structures."""

# python 2/3 compatiblity
from __future__ import division, print_function, absolute_import

# global imports
from lxml import etree

__all__ = ['is_true', 'MachineryComposition', 'SpeciesReference',
           'ListOfReactants', 'ListOfProducts', 'TargetValue']


def is_true(attribute):
    """
    Determine whether an XML attribute evaluates as the logical true value.

    Parameters
    ----------
    attribute : str
        attribute value.

    Returns
    -------
    out : bool
        True if attribute value is considered to represent the 'true' boolean
        state, False otherwise.

    """
    return attribute.lower() == 'true' or attribute == '1'


def get_unique_child(parent, child_name, strict=True):
    """
    Get children node matching given name.

    Parameters
    ----------
    parent : xml element
        Parent node.
    child_name : str
        Name of the node to look for.
    strict : bool, optional
        Flag telling whether there should be strictly one
        child. If set to True, an error is raised if no child was found.

    Returns
    -------
    out : xml element or None
        Child if exactly one node is found. If no child is found,
        None is returned if strict was set to False.

    Raises
    ------
    UserWarning
        When no child is found and strict is set to True or
        when more than one node is found.

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
    List of objects of same type.

    Attributes
    ----------
    list_element : class
        type of elements stored by list.

    """

    list_element = None

    def __init__(self):
        """Build default object."""
        self._elements = []
        self._elements_by_id = {}

    def __getitem__(self, i):
        """Return item with given index."""
        return self._elements[i]

    def __iter__(self):
        """Return iterator over list."""
        return iter(self._elements)

    def __len__(self):
        """Return length of list."""
        return len(self._elements)

    def append(self, element):
        """Append element to list."""
        assert isinstance(element, self.list_element)
        self._elements.append(element)
        try:
            self._elements_by_id[element.id] = element
        except AttributeError:
            pass

    def remove(self, element):
        """Remove element from list."""
        self._elements.remove(element)
        try:
            del self._elements_by_id[element.id]
        except AttributeError:
            pass

    def get_by_id(self, identifier):
        """Return element with given identifier or None if not found."""
        return self._elements_by_id.get(identifier)

    def is_empty(self):
        """Return whether list is empty."""
        return len(self._elements) == 0

    def to_xml_node(self):
        """Convert to xml node."""
        result = etree.Element(self.tag)
        result.extend([e.to_xml_node() for e in self._elements])
        return result

    @classmethod
    def from_xml_node(cls, node):
        """Build object from xml node."""
        result = cls()
        for n in node.iterfind(cls.list_element.tag):
            result.append(cls.list_element.from_xml_node(n))
        return result


class MachineryComposition(object):
    """
    Machinery composition information.

    Attributes
    ----------
    reactants : ListOfReactants
        List of reactants used to assemble machinery.
    products :ListOfProducts
        List of byproducts generated while assembling machinery.

    """

    tag = 'machineryComposition'

    def __init__(self):
        """Build default object."""
        self.reactants = ListOfReactants()
        self.products = ListOfProducts()

    def is_empty(self):
        """Return whether composition is fully empty."""
        return self.reactants.is_empty() and self.products.is_empty()

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


class SpeciesReference(object):
    """
    Reference to a species and its stoichiometry.

    Attributes
    ----------
    species : str
        Identifier of species.
    stoichiometry : float
        Stoichiometry of species.

    """

    tag = 'speciesReference'

    def __init__(self, species, stoichiometry):
        """
        Constructor.

        Parameters
        ----------
        species : str
            Identifier of species.
        stoichiometry : float
            Stoichiometry of species.

        """
        self.species = species
        self.stoichiometry = stoichiometry

    def to_xml_node(self):
        """Convert to xml node."""
        result = etree.Element(self.tag)
        result.set('species', self.species)
        result.set('stoichiometry', str(self.stoichiometry))
        return result

    @classmethod
    def from_xml_node(cls, node):
        """Build object from xml node."""
        return cls(node.get('species'), float(node.get('stoichiometry')))


class ListOfReactants(ListOf):
    """List of SpeciesReference representing reactants."""

    tag = 'listOfReactants'
    list_element = SpeciesReference


class ListOfProducts(ListOf):
    """List of SpeciesReference representing products."""

    tag = 'listOfProducts'
    list_element = SpeciesReference


class TargetValue(object):
    """
    Target specifying value, lower bound and/or upper bound of a variable.

    Attributes
    ----------
    value : float or None
        exact target value (None if no exact value to match).
    lower_bound : float or None
        lower bound on target value (None if no lower bound).
    upper_bound : float or None
        upper bound on target value (None if no upper bound).

    """

    tag = 'targetValue'

    def __init__(self):
        """Build default object."""
        self.value = None
        self.lower_bound = None
        self.upper_bound = None

    def is_empty(self):
        """Return whether all attributes are unspecified."""
        return self.value is None and self.lower_bound is None \
            and self.upper_bound is None

    def to_xml_node(self):
        """Convert to xml node."""
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
        """Build object from xml node."""
        result = cls()
        result._init_from_xml_node(node)
        return result

    def _init_from_xml_node(self, node):
        """Match attributes with given node."""
        self.value = node.get('value')
        self.lower_bound = node.get('lowerBound')
        self.upper_bound = node.get('upperBound')
