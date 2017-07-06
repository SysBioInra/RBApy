"""
Module defining process-specific classes used for RBA XML structures.
"""

# python 2/3 compatiblity
from __future__ import division, print_function

# global imports
from lxml import etree

# local imports
from rba.xml.common import (get_unique_child, ListOf,
                            TargetValue, MachineryComposition,
                            ListOfReactants, ListOfProducts)

__all__ = ['RbaProcesses', 'Process', 'ListOfProcesses',
           'Machinery', 'Capacity', 'Operations', 'Operation',
           'ListOfProductions', 'ListOfDegradations',
           'Targets', 'TargetSpecies', 'TargetReaction',
           'ListOfConcentrations', 'ListOfProductionFluxes',
           'ListOfDegradationFluxes', 'ListOfReactionFluxes',
           'ComponentMap', 'ListOfComponentMaps',
           'ConstantCost', 'Cost', 'ListOfCosts']


class RbaProcesses(object):
    def __init__(self):
        self.processes = ListOfProcesses()
        self.component_maps = ListOfComponentMaps()

    def write(self, output_stream, doc_name = 'RBAProcesses'):
        root = etree.Element(doc_name)
        root.extend([self.processes.to_xml_node(),
                     self.component_maps.to_xml_node()])
        etree.ElementTree(root).write(output_stream, pretty_print=True)

    @classmethod
    def from_file(cls, input_stream):
        node = etree.ElementTree(file=input_stream).getroot()
        result = cls()
        n = get_unique_child(node, ListOfProcesses.tag)
        result.processes = ListOfProcesses.from_xml_node(n)
        n = get_unique_child(node, ListOfComponentMaps.tag)
        result.component_maps = ListOfComponentMaps.from_xml_node(n)
        return result

class Process(object):
    tag = 'process'
    
    def __init__(self, id_, name):
        self.id = id_
        self.name = name
        self.machinery = Machinery()
        self.operations = Operations()
        self.targets = Targets()

    def to_xml_node(self):
        process = etree.Element(self.tag)
        process.set('id', self.id)
        process.set('name', self.name)
        # optional subelements
        opt = [self.machinery, self.operations, self.targets]
        process.extend([e.to_xml_node() for e in opt if not(e.is_empty())])
        return process

    @classmethod
    def from_xml_node(cls, node):
        result = cls(node.get('id'), node.get('name'))
        n = get_unique_child(node, Machinery.tag, False)
        if n is not None: result.machinery = Machinery.from_xml_node(n)
        n = get_unique_child(node, Operations.tag, False)
        if n is not None: result.operations = Operations.from_xml_node(n)
        n = get_unique_child(node, Targets.tag, False)
        if n is not None: result.targets = Targets.from_xml_node(n)
        return result

class ListOfProcesses(ListOf):
    tag = 'listOfProcesses'
    list_element = Process

class Machinery(object):
    tag = 'machinery'
    
    def __init__(self):
        self.machinery_composition = MachineryComposition()
        self.capacity = Capacity()

    def is_empty(self):
        return self.machinery_composition.is_empty() \
            and self.capacity.is_empty()

    def to_xml_node(self):
        result = etree.Element(self.tag)
        result.extend([self.machinery_composition.to_xml_node(),
                       self.capacity.to_xml_node()])
        return result

    @classmethod
    def from_xml_node(cls, node):
        result = cls()
        machinery_node = get_unique_child(node, MachineryComposition.tag)
        result.machinery_composition \
            = MachineryComposition.from_xml_node(machinery_node)
        capacity_node = get_unique_child(node, Capacity.tag)
        result.capacity = Capacity.from_xml_node(capacity_node)
        return result

class Capacity(TargetValue):
    tag = 'capacity'

class Operations(object):
    tag = 'operations'
    
    def __init__(self):
        self.productions = ListOfProductions()
        self.degradations = ListOfDegradations()

    def is_empty(self):
        return self.productions.is_empty() and self.degradations.is_empty()

    def to_xml_node(self):
        result = etree.Element(self.tag)
        # optional subelements
        opt = [self.productions, self.degradations]
        result.extend([e.to_xml_node() for e in opt if not(e.is_empty())])
        return result
        
    @classmethod
    def from_xml_node(cls, node):
        result = cls()
        n = get_unique_child(node, ListOfProductions.tag, False)
        if n is not None:
            result.productions = ListOfProductions.from_xml_node(n)
        n = get_unique_child(node, ListOfDegradations.tag, False)
        if n is not None:
            result.degradations = ListOfDegradations.from_xml_node(n)
        return result

class Operation(object):
    tag = 'operation'
    
    def __init__(self, map_, set_):
        self.component_map = map_
        self.set = set_

    def to_xml_node(self):
        result = etree.Element(self.tag)
        result.set('componentMap', self.component_map)
        result.set('set', self.set)
        return result

    @classmethod
    def from_xml_node(cls, node):
        result = cls('', '')
        result.component_map = node.get('componentMap')
        result.set = node.get('set')
        return result

class ListOfProductions(ListOf):
    tag = 'listOfProductions'
    list_element = Operation

class ListOfDegradations(ListOf):
    tag = 'listOfDegradations'
    list_element = Operation

class Targets(object):
    tag = 'targets'
    
    def __init__(self):
        self.concentrations = ListOfConcentrations()
        self.production_fluxes = ListOfProductionFluxes()
        self.degradation_fluxes = ListOfDegradationFluxes()
        self.reaction_fluxes = ListOfReactionFluxes()

    def is_empty(self):
        return all((l.is_empty() for l in [self.concentrations, self.production_fluxes, self.degradation_fluxes, self.reaction_fluxes]))

    def to_xml_node(self):
        result = etree.Element(self.tag)
        # optional subelements
        opt = [self.concentrations, self.production_fluxes,
               self.degradation_fluxes, self.reaction_fluxes]
        result.extend([e.to_xml_node() for e in opt if not(e.is_empty())])
        return result

    @classmethod
    def from_xml_node(cls, node):
        result = cls()
        n = get_unique_child(node, result.concentrations.tag, False)
        if n is not None:
            result.concentrations = result.concentrations.from_xml_node(n)
        n = get_unique_child(node, result.production_fluxes.tag, False)
        if n is not None:
            result.production_fluxes = result.production_fluxes.from_xml_node(n)
        n = get_unique_child(node, result.degradation_fluxes.tag, False)
        if n is not None:
            result.degradation_fluxes = result.degradation_fluxes.from_xml_node(n)
        n = get_unique_child(node, result.reaction_fluxes.tag, False)
        if n is not None:
            result.reaction_fluxes = result.reaction_fluxes.from_xml_node(n)
        return result
    
class TargetSpecies(TargetValue):
    tag = 'targetSpecies'
    
    def __init__(self, species):
        super(TargetSpecies, self).__init__()
        self.species = species

    def to_xml_node(self):
        result = super(TargetSpecies, self).to_xml_node()
        result.set('species', self.species)
        return result

    @classmethod
    def from_xml_node(cls, node):
        result = cls(node.get('species'))
        result._init_from_xml_node(node)
        return result

class TargetReaction(TargetValue):
    tag = 'targetReaction'
    
    def __init__(self, reaction):
        super(TargetReaction, self).__init__()
        self.reaction = reaction

    def to_xml_node(self):
        result = super(TargetReaction, self).to_xml_node()
        result.set('reaction', self.reaction)
        return result

    @classmethod
    def from_xml_node(cls, node):
        result = cls(node.get('reaction'))
        result._init_from_xml_node(node)
        return result

class ListOfConcentrations(ListOf):
    tag = 'listOfConcentrations'
    list_element = TargetSpecies

class ListOfProductionFluxes(ListOf):
    tag = 'listOfProductionFluxes'
    list_element = TargetSpecies

class ListOfDegradationFluxes(ListOf):
    tag = 'listOfDegradationFluxes'
    list_element = TargetSpecies

class ListOfReactionFluxes(ListOf):
    tag = 'listOfReactionFluxes'
    list_element = TargetReaction

class ComponentMap(object):
    tag = 'componentMap'
    
    def __init__(self, id_):
        self.id = id_
        self.constant_cost = ConstantCost()
        self.costs = ListOfCosts()
        
    def to_xml_node(self):
        result = etree.Element(self.tag)
        result.set('id', self.id)
        result.extend([self.constant_cost.to_xml_node(),
                       self.costs.to_xml_node()])
        return result

    @classmethod
    def from_xml_node(cls, node):
        result = cls(node.get('id'))
        n = get_unique_child(node, ConstantCost.tag, False)
        if n is not None: result.constant_cost = ConstantCost.from_xml_node(n)
        n = get_unique_child(node, ListOfCosts.tag)
        result.costs = ListOfCosts.from_xml_node(n)
        return result

class ListOfComponentMaps(ListOf):
    tag = 'listOfComponentMaps'
    list_element = ComponentMap

class ConstantCost(object):
    tag = 'constantCost'
    
    def __init__(self):
        self.reactants = ListOfReactants()
        self.products = ListOfProducts()

    def to_xml_node(self):
        result = etree.Element(self.tag)
        # optional subelements
        opt = [self.reactants, self.products]
        result.extend([e.to_xml_node() for e in opt if not(e.is_empty())])
        return result

    @classmethod
    def from_xml_node(cls, node):
        result = cls()
        n = get_unique_child(node, ListOfReactants.tag, False)
        if n is not None: result.reactants = ListOfReactants.from_xml_node(n)
        n = get_unique_child(node, ListOfProducts.tag, False)
        if n is not None: result.products = ListOfProducts.from_xml_node(n)
        return result

class Cost(object):
    tag = 'cost'
    
    def __init__(self, component, processing_cost=0):
        self.component = component
        self.processing_cost = processing_cost
        self.reactants = ListOfReactants()
        self.products = ListOfProducts()

    def to_xml_node(self):
        result = etree.Element(self.tag)
        result.set('component', self.component)
        result.set('processingCost', str(self.processing_cost))
        # optional subelements
        opt = [self.reactants, self.products]
        result.extend([e.to_xml_node() for e in opt if not(e.is_empty())])
        return result

    @classmethod
    def from_xml_node(cls, node):
        try:
            proc_cost = float(node.get('processingCost'))
        except TypeError:
            proc_cost = 0
        result = cls(node.get('component'), proc_cost)
        n = get_unique_child(node, ListOfReactants.tag, False)
        if n is not None: result.reactants = ListOfReactants.from_xml_node(n)
        n = get_unique_child(node, ListOfProducts.tag, False)
        if n is not None: result.products = ListOfProducts.from_xml_node(n)
        return result

class ListOfCosts(ListOf):
    tag = 'listOfCosts'
    list_element = Cost
