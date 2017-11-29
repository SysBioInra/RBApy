"""Module defining target-specific classes used for RBA XML structures."""

# python 2/3 compatiblity
from __future__ import division, print_function, absolute_import

# global imports
from lxml import etree

# local imports
from rba.xml.common import get_unique_child, ListOf, TargetValue

__all__ = ['RbaTargets', 'TargetGroup', 'ListOfTargetGroups',
           'TargetSpecies', 'TargetReaction',
           'ListOfConcentrations', 'ListOfProductionFluxes',
           'ListOfDegradationFluxes', 'ListOfReactionFluxes']


class RbaTargets(object):
    """
    Target-related structures.

    Attributes
    ----------
    targets : Targets
        Metabolic or reaction fluxes maintained/induced by process.

    """

    def __init__(self):
        """Constructor."""
        self.target_groups = ListOfTargetGroups()

    def write(self, output_stream, doc_name='RBATargets'):
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
        root.extend([self.target_groups.to_xml_node()])
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
        n = get_unique_child(node, ListOfTargetGroups.tag)
        result.target_groups = ListOfTargetGroups.from_xml_node(n)
        return result


class TargetGroup(object):
    """
    Group of target fluxes.

    Attributes
    ----------
    id : str
        Identifier.
    concentrations : ListOfConcentrations
        Concentration of metabolites to keep at given level.
    production_fluxes : ListOfProductionFluxes
        Metabolite/macromolecule production fluxes.
    degradation_fluxes : ListOfDegradationFluxes
        Macromolecule degradation/metabolite consumption fluxes.
    reaction_fluxes : ListOfReactionFluxes
        Target fluxes for metabolic reactions.

    """

    tag = 'targetGroup'

    def __init__(self, id_):
        """
        Constructor.

        Parameters
        ----------
        id_ : str
            Identifier.

        """
        self.id = id_
        self.concentrations = ListOfConcentrations()
        self.production_fluxes = ListOfProductionFluxes()
        self.degradation_fluxes = ListOfDegradationFluxes()
        self.reaction_fluxes = ListOfReactionFluxes()

    def is_empty(self):
        """Return whether any target was specified."""
        return all((l.is_empty()
                    for l in [self.concentrations, self.production_fluxes,
                              self.degradation_fluxes, self.reaction_fluxes]))

    def to_xml_node(self):
        """Convert to xml node."""
        result = etree.Element(self.tag)
        result.set('id', self.id)
        # optional subelements
        opt = [self.concentrations, self.production_fluxes,
               self.degradation_fluxes, self.reaction_fluxes]
        result.extend([e.to_xml_node() for e in opt if not e.is_empty()])
        return result

    @classmethod
    def from_xml_node(cls, node):
        """Build object from xml node."""
        result = cls(node.get('id'))
        n = get_unique_child(node, result.concentrations.tag, False)
        if n is not None:
            result.concentrations = result.concentrations.from_xml_node(n)
        n = get_unique_child(node, result.production_fluxes.tag, False)
        if n is not None:
            result.production_fluxes \
                = result.production_fluxes.from_xml_node(n)
        n = get_unique_child(node, result.degradation_fluxes.tag, False)
        if n is not None:
            result.degradation_fluxes \
                = result.degradation_fluxes.from_xml_node(n)
        n = get_unique_child(node, result.reaction_fluxes.tag, False)
        if n is not None:
            result.reaction_fluxes = result.reaction_fluxes.from_xml_node(n)
        return result


class ListOfTargetGroups(ListOf):
    """Groups used to define target fluxes."""

    tag = 'listOfTargetGroups'
    list_element = TargetGroup


class TargetSpecies(TargetValue):
    """
    Association of a target value with a chemical species.

    Attributes
    ----------
    species : str
        Identifier of species targetted.

    """

    tag = 'targetSpecies'

    def __init__(self, species):
        """
        Constructor.

        Parameters
        ----------
        species : str
            Identifier of species.

        """
        super(TargetSpecies, self).__init__()
        self.species = species

    def to_xml_node(self):
        """Convert to xml node."""
        result = super(TargetSpecies, self).to_xml_node()
        result.set('species', self.species)
        return result

    @classmethod
    def from_xml_node(cls, node):
        """Build object from xml node."""
        result = cls(node.get('species'))
        result._init_from_xml_node(node)
        return result


class TargetReaction(TargetValue):
    """
    Association of a target value with a chemical reaction.

    Attributes
    ----------
    reaction : str
        Identifier of reaction targetted.

    """

    tag = 'targetReaction'

    def __init__(self, reaction):
        """
        Constructor.

        Parameters
        ----------
        reaction : str
            Reaction identifier.

        """
        super(TargetReaction, self).__init__()
        self.reaction = reaction

    def to_xml_node(self):
        """Convert to xml node."""
        result = super(TargetReaction, self).to_xml_node()
        result.set('reaction', self.reaction)
        return result

    @classmethod
    def from_xml_node(cls, node):
        """Build object from xml node."""
        result = cls(node.get('reaction'))
        result._init_from_xml_node(node)
        return result


class ListOfConcentrations(ListOf):
    """Metabolic concentrations to maintain at a given level."""

    tag = 'listOfConcentrations'
    list_element = TargetSpecies


class ListOfProductionFluxes(ListOf):
    """Production fluxes of metabolites/macromolecules."""

    tag = 'listOfProductionFluxes'
    list_element = TargetSpecies


class ListOfDegradationFluxes(ListOf):
    """Fluxes of metabolite cosumption/macromolecule degradation."""

    tag = 'listOfDegradationFluxes'
    list_element = TargetSpecies


class ListOfReactionFluxes(ListOf):
    """Fluxes of metabolic reactions to maintain at given level."""

    tag = 'listOfReactionFluxes'
    list_element = TargetReaction


class ComponentMap(object):
    """
    Cost of [un]polymerization of components of a macromolecule.

    Attributes
    ----------
    id : str
        Identifier of map.
    constant_cost : ConstantCost
        Cost associated with synthesis/degradation that do not depend
        on the number of components (e.g. initialization).
    cost : ListOfCosts
        Cost associated with [un]polymerization of components.

    """

    tag = 'componentMap'

    def __init__(self, id_):
        """
        Constructor.

        Parameters
        ----------
        id_ : str
            Identifier of map.

        """
        self.id = id_
        self.constant_cost = ConstantCost()
        self.costs = ListOfCosts()

    def to_xml_node(self):
        """Convert to xml node."""
        result = etree.Element(self.tag)
        result.set('id', self.id)
        result.extend([self.constant_cost.to_xml_node(),
                       self.costs.to_xml_node()])
        return result

    @classmethod
    def from_xml_node(cls, node):
        """Build object from xml node."""
        result = cls(node.get('id'))
        n = get_unique_child(node, ConstantCost.tag, False)
        if n is not None:
            result.constant_cost = ConstantCost.from_xml_node(n)
        n = get_unique_child(node, ListOfCosts.tag)
        result.costs = ListOfCosts.from_xml_node(n)
        return result


class ListOfComponentMaps(ListOf):
    """List of ComponentMap elements."""

    tag = 'listOfComponentMaps'
    list_element = ComponentMap


class ConstantCost(object):
    """
    Metabolites consumed/produced independent of component number/nature.

    Attributes
    ----------
    reactants : ListOfReactants
        Metabolites consumed.
    products : ListOfProducts
        Metabolites generated as byproducts.

    """

    tag = 'constantCost'

    def __init__(self):
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


class Cost(object):
    """
    Cost of [un]polymerization of given component.

    Attributes
    ----------
    component : str
        Identifier of component.
    processing_cost : float
        Capacity cost for process machinery to [un]polymerize component.
    reactants : ListOfReactants
        Metabolites consumed.
    products : ListOfProducts
        Metabolites generated as byproducts.

    """

    tag = 'cost'

    def __init__(self, component, processing_cost=0.0):
        """
        Constructor.

        Parameters
        ----------
        component : str
            Identifier of component.
        processing_cost : float, optional
            Capacity cost for process machinery to [un]polymerize component.

        """
        self.component = component
        self.processing_cost = processing_cost
        self.reactants = ListOfReactants()
        self.products = ListOfProducts()

    def to_xml_node(self):
        """Convert to xml node."""
        result = etree.Element(self.tag)
        result.set('component', self.component)
        result.set('processingCost', str(self.processing_cost))
        # optional subelements
        opt = [self.reactants, self.products]
        result.extend([e.to_xml_node() for e in opt if not e.is_empty()])
        return result

    @classmethod
    def from_xml_node(cls, node):
        """Build object from xml node."""
        try:
            proc_cost = float(node.get('processingCost'))
        except TypeError:
            proc_cost = 0.0
        result = cls(node.get('component'), proc_cost)
        n = get_unique_child(node, ListOfReactants.tag, False)
        if n is not None:
            result.reactants = ListOfReactants.from_xml_node(n)
        n = get_unique_child(node, ListOfProducts.tag, False)
        if n is not None:
            result.products = ListOfProducts.from_xml_node(n)
        return result


class ListOfCosts(ListOf):
    """List of Cost elements."""

    tag = 'listOfCosts'
    list_element = Cost
