
from lxml import etree
#import xml.etree.cElementTree as etree

################################################################################
#                                  COMMON                                      #
################################################################################

def is_true(attribute):
    """
    Determine whether an XML attribute evaluates as the logical true value.
    """
    return attribute.lower() == 'true' or attribute == '1'

def get_unique_child(parent, child_name, strict = True):
    """
    Function looking for a children node matching a given name.
    :param parent: Parent node.
    :param child_name: Name of the node to look for.
    :param strict: Boolean flag telling whether there should be strictly one
     child. If set to True, an error is raised if no child was found.
    :return: Element if exactly one node is found. If no child is found,
     a UserWarning exception is raised if strict is set to True, 
     otherwise None is returned.
     A UserWarning exception is raised if more than one node is found.
    """
    children = parent.findall(child_name)
    if len(children) == 1:
        return children[0]
    elif (len(children) == 0) and not(strict):
        return None
    else:
        if strict:
            raise UserWarning('A ' + parent.tag + ' should have exactly '
                              'one ' + child_name + ' node child.')
        else:
            raise UserWarning('A ' + parent.tag + ' should not have more '
                              'than one ' + child_name + ' node children.')

class MachineryComposition(object):
    tag = 'machineryComposition'
    
    def __init__(self):
        self.reactants = ListOfReactants()
        self.products = ListOfProducts()

    def is_empty(self):
        return self.reactants.is_empty() and self.products.is_empty()

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
        if n is not None:
            result.reactants = ListOfReactants.from_xml_node(n)
        n = get_unique_child(node, ListOfProducts.tag, False)
        if n is not None:
            result.products = ListOfProducts.from_xml_node(n)
        return result
    
class ListOf(object):
    def __init__(self):
        self._elements = []

    def __getitem__(self, i):
        return self._elements[i]

    def __len__(self):
        return len(self._elements)

    def append(self, element):
        assert(isinstance(element, self.list_element))
        self._elements.append(element)

    def is_empty(self):
        return len(self._elements) == 0

    def to_xml_node(self):
        result = etree.Element(self.tag)
        result.extend([e.to_xml_node() for e in self._elements])
        return result

    @classmethod
    def from_xml_node(cls, node):
        result = cls()
        for n in node.iterfind(cls.list_element.tag):
            result.append(cls.list_element.from_xml_node(n))
        return result

class SpeciesReference(object):
    tag = 'speciesReference'
    
    def __init__(self, species, stoichiometry):
        self.species = species
        self.stoichiometry = stoichiometry

    def to_xml_node(self):
        result = etree.Element(self.tag)
        result.set('species', self.species)
        result.set('stoichiometry', str(self.stoichiometry))
        return result

    @classmethod
    def from_xml_node(cls, node):
        return cls(node.get('species'), float(node.get('stoichiometry')))

class ListOfReactants(ListOf):
    tag = 'listOfReactants'
    list_element = SpeciesReference
        
class ListOfProducts(ListOf):
    tag = 'listOfProducts'
    list_element = SpeciesReference

class Parameter(object):
    tag = 'parameter'
    
    def __init__(self, id_, value):
        self.id = id_
        self.value = value

    def to_xml_node(self):
        result = etree.Element(self.tag)
        result.set('id', self.id)
        result.set('value', str(self.value))
        return result

    @classmethod
    def from_xml_node(cls, node):
        return cls(node.get('id'), float(node.get('value')))

class ListOfParameters(ListOf):
    tag = 'listOfParameters'
    list_element = Parameter

class Function(object):
    tag = 'function'
    
    def __init__(self, id_, type_, parameters = {}, variable = ''):
        if id_ is not None: self.id = id_
        else: self.id = ''
        self.type = type_
        self.variable = variable
        self.parameters = ListOfParameters()
        for key, value in parameters.items():
            self.parameters.append(Parameter(key, value))

    def to_xml_node(self):
        result = etree.Element(self.tag)
        result.set('id', self.id)
        result.set('type', self.type)
        if self.variable: result.set('variable', self.variable)
        if not(self.parameters.is_empty()):
            result.extend([self.parameters.to_xml_node()])
        return result

    @classmethod
    def from_xml_node(cls, node):
        result = cls(node.get('id'), node.get('type'), {}, node.get('variable'))
        n = get_unique_child(node, 'listOfParameters', False)
        if n is not None:
            result.parameters = ListOfParameters.from_xml_node(n)
        return result

class ListOfFunctions(ListOf):
    tag = 'listOfFunctions'
    list_element = Function

class TargetValue(object):
    tag = 'targetValue'
    
    def __init__(self):
        self.value = None
        self.lower_bound = None
        self.upper_bound = None

    def is_empty(self):
        return self.value is None and self.lower_bound is None \
            and self.upper_bound is None

    def to_xml_node(self):
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
        result = cls()
        result._init_from_xml_node(node)
        return result

    def _init_from_xml_node(self, node):
        self.value = node.get('value')
        self.lower_bound = node.get('lowerBound')
        self.upper_bound = node.get('upperBound')

################################################################################
#                                  METABOLISM                                  #
################################################################################

class RbaMetabolism(object):
    def __init__(self):
        self.compartments = ListOfCompartments()
        self.species = ListOfSpecies()
        self.reactions = ListOfReactions()

    def write(self, output_stream, doc_name = 'RBAMetabolism'):
        root = etree.Element(doc_name)
        root.extend([self.compartments.to_xml_node(),
                     self.species.to_xml_node(), self.reactions.to_xml_node()])
        etree.ElementTree(root).write(output_stream, pretty_print=True)

    @classmethod
    def from_file(cls, input_stream):
        node = etree.ElementTree(file=input_stream).getroot()
        result = cls()
        n = get_unique_child(node, ListOfCompartments.tag)
        result.compartments = ListOfCompartments.from_xml_node(n)
        n = get_unique_child(node, ListOfSpecies.tag)
        result.species = ListOfSpecies.from_xml_node(n)
        n = get_unique_child(node, ListOfReactions.tag)
        result.reactions = ListOfReactions.from_xml_node(n)
        return result

class Compartment(object):
    tag = 'compartment'
    
    def __init__(self, id_):
        self.id = id_

    def to_xml_node(self):
        result = etree.Element(self.tag)
        result.set('id', self.id)
        return result

    @classmethod
    def from_xml_node(cls, node):
        return cls(node.get('id'))
    
class ListOfCompartments(ListOf):
    tag = 'listOfCompartments'
    list_element = Compartment

class Species(object):
    tag = 'species'
    
    def __init__(self, id_, boundary_condition):
        self.id = id_
        self.boundary_condition = boundary_condition

    def to_xml_node(self):
        result = etree.Element(self.tag)
        result.set('id', self.id)
        result.set('boundaryCondition', str(self.boundary_condition).lower())
        return result

    @classmethod
    def from_xml_node(cls, node):
        return cls(node.get('id'), is_true(node.get('boundaryCondition')))

class ListOfSpecies(ListOf):
    tag = 'listOfSpecies'
    list_element = Species

class Reaction(object):
    tag = 'reaction'
    
    def __init__(self, id_, reversible):
        self.id = id_
        self.reversible = reversible
        self.reactants = ListOfReactants()
        self.products = ListOfProducts()

    def to_xml_node(self):
        result = etree.Element(self.tag)
        result.set('id', self.id)
        result.set('reversible', str(self.reversible).lower())
        result.extend([self.reactants.to_xml_node(),
                       self.products.to_xml_node()])
        return result

    @classmethod
    def from_xml_node(cls, node):
        result = cls(node.get('id'), is_true(node.get('reversible')))
        n = get_unique_child(node, ListOfReactants.tag, False)
        if n is not None: result.reactants = ListOfReactants.from_xml_node(n)
        n = get_unique_child(node, ListOfProducts.tag, False)
        if n is not None: result.products = ListOfProducts.from_xml_node(n)
        return result

class ListOfReactions(ListOf):
    tag = 'listOfReactions'
    list_element = Reaction

################################################################################
#                              PARAMETERS                                      #
################################################################################

class RbaParameters(object):
    def __init__(self):
        self.target_densities = ListOfTargetDensities()
        self.functions = ListOfFunctions()
        self.aggregates = ListOfAggregates()

    def write(self, output_stream, doc_name = 'RBAParameters'):
        root = etree.Element(doc_name)
        root.extend([self.target_densities.to_xml_node(),
                     self.functions.to_xml_node(),
                     self.aggregates.to_xml_node()])
        etree.ElementTree(root).write(output_stream, pretty_print=True)

    @classmethod
    def from_file(cls, input_stream):
        node = etree.ElementTree(file=input_stream).getroot()
        result = cls()
        n = get_unique_child(node, ListOfTargetDensities.tag)
        result.target_densities = ListOfTargetDensities.from_xml_node(n)
        n = get_unique_child(node, ListOfFunctions.tag)
        result.functions = ListOfFunctions.from_xml_node(n)
        n = get_unique_child(node, ListOfAggregates.tag)
        result.aggregates = ListOfAggregates.from_xml_node(n)
        return result
    
class TargetDensity(TargetValue):
    tag = 'targetDensity'
    
    def __init__(self, compartment):
        super(TargetDensity, self).__init__()
        self.compartment = compartment

    def to_xml_node(self):
        result = super(TargetDensity,self).to_xml_node()
        result.set('compartment', self.compartment)
        return result
    
    @classmethod
    def from_xml_node(cls, node):
        result = cls(node.get('compartment'))
        result._init_from_xml_node(node)
        return result

class ListOfTargetDensities(ListOf):
    tag = 'listOfTargetDensities'
    list_element = TargetDensity

class FunctionReference(object):
    tag = 'functionReference'

    def __init__(self, function):
        self.function = function

    def to_xml_node(self):
        result = etree.Element(self.tag)
        result.set('function', self.function)
        return result

    @classmethod
    def from_xml_node(cls, node):
        return cls(node.get('function'))

class ListOfFunctionReferences(ListOf):
    tag = 'listOfFunctionReferences'
    list_element = FunctionReference
    
class Aggregate(object):
    tag = 'aggregate'

    def __init__(self, id_, type_):
        if id_ is not None: self.id = id_
        else: self.id = ''
        self.type = type_
        self.function_references = ListOfFunctionReferences()

    def to_xml_node(self):
        result = etree.Element(self.tag)
        result.set('id', self.id)
        result.set('type', self.type)
        result.extend([self.function_references.to_xml_node()])
        return result

    @classmethod
    def from_xml_node(cls, node):
        result = cls(node.get('id'), node.get('type'))
        n = get_unique_child(node, ListOfFunctionReferences.tag)
        result.function_references = ListOfFunctionReferences.from_xml_node(n)
        return result

class ListOfAggregates(ListOf):
    tag = 'listOfAggregates'
    list_element = Aggregate

################################################################################
#                               MACROMOLECULES                                 #
################################################################################

class RbaMacromolecules(object):
    def __init__(self):
        self.components = ListOfComponents()
        self.macromolecules = ListOfMacromolecules()

    def write(self, output_stream, doc_name = 'RBAMacromolecules'):
        root = etree.Element(doc_name)
        root.extend([self.components.to_xml_node(),
                     self.macromolecules.to_xml_node()])
        etree.ElementTree(root).write(output_stream, pretty_print=True)

    @classmethod
    def from_file(cls, input_stream):
        node = etree.ElementTree(file=input_stream).getroot()
        result = cls()
        n = get_unique_child(node, ListOfComponents.tag)
        result.components = ListOfComponents.from_xml_node(n)
        n = get_unique_child(node, ListOfMacromolecules.tag)
        result.macromolecules = ListOfMacromolecules.from_xml_node(n)
        return result

class Component(object):
    tag = 'component'
    
    def __init__(self, id_, name, type_, weight):
        self.id = id_
        self.name = name
        self.type = type_
        self.weight = weight

    def to_xml_node(self):
        result = etree.Element(self.tag)
        result.set('id', self.id)
        result.set('name', self.name)
        result.set('type', self.type)
        result.set('weight', str(self.weight))
        return result

    @classmethod
    def from_xml_node(cls, node):
        return cls(node.get('id'), node.get('name'),
                   node.get('type'), node.get('weight'))

class ListOfComponents(ListOf):
    tag = 'listOfComponents'
    list_element = Component

class Macromolecule(object):
    tag = 'macromolecule'
    
    def __init__(self, id_, compartment, composition = {}):
        self.id = id_
        self.compartment = compartment
        self.composition = Composition()
        for comp, sto in composition.items():
            self.composition.append(ComponentReference(comp, sto))

    def to_xml_node(self):
        result = etree.Element(self.tag)
        result.set('id', self.id)
        result.set('compartment', self.compartment)
        if not(self.composition.is_empty()):
            result.extend([self.composition.to_xml_node()])
        return result

    @classmethod
    def from_xml_node(cls, node):
        result = cls(node.get('id'), node.get('compartment'))
        comp_node = get_unique_child(node, Composition.tag)
        result.composition = Composition.from_xml_node(comp_node)
        return result

class ListOfMacromolecules(ListOf):
    tag = 'listOfMacromolecules'
    list_element = Macromolecule

class ComponentReference(object):
    tag = 'componentReference'
    
    def __init__(self, component, stoichiometry):
        self.component = component
        self.stoichiometry = stoichiometry

    def to_xml_node(self):
        result = etree.Element(self.tag)
        result.set('component', self.component)
        result.set('stoichiometry', str(self.stoichiometry))
        return result

    @classmethod
    def from_xml_node(cls, node):
        return cls(node.get('component'), float(node.get('stoichiometry')))

class Composition(ListOf):
    tag = 'composition'
    list_element = ComponentReference

################################################################################
#                                  PROCESSES                                   #
################################################################################

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
    
    def __init__(self, component, processing_cost = 0):
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

################################################################################
#                                   ENZYMES                                    #
################################################################################

class RbaEnzymes(object):
    def __init__(self):
        self.efficiency_functions = ListOfEfficiencyFunctions()
        self.enzymes = ListOfEnzymes()

    def write(self, output_stream, doc_name = 'RBAEnzymes'):
        root = etree.Element(doc_name)
        root.extend([self.efficiency_functions.to_xml_node(),
                     self.enzymes.to_xml_node()])
        etree.ElementTree(root).write(output_stream, pretty_print=True)

    @classmethod
    def from_file(cls, input_stream):
        node = etree.ElementTree(file=input_stream).getroot()
        result = cls()
        n = get_unique_child(node, ListOfEfficiencyFunctions.tag)
        result.efficiency_functions = ListOfEfficiencyFunctions.from_xml_node(n)
        n = get_unique_child(node, ListOfEnzymes.tag)
        result.enzymes = ListOfEnzymes.from_xml_node(n)
        return result

class ListOfEfficiencyFunctions(ListOf):
    tag = 'listOfEfficiencyFunctions'
    list_element = Function

class Enzyme(object):
    tag = 'enzyme'
    
    def __init__(self, reaction, zero_cost = False):
        if zero_cost is None: zero_cost = False
        self.id = reaction
        self.zero_cost = zero_cost
        self.machinery_composition = MachineryComposition()
        self.enzymatic_activity = EnzymaticActivity(reaction)

    def to_xml_node(self):
        result = etree.Element(self.tag)
        result.set('id', self.id)
        result.set('zeroCost', str(self.zero_cost).lower())
        if not(self.machinery_composition.is_empty()):
            result.extend([self.machinery_composition.to_xml_node()])
        result.extend([self.enzymatic_activity.to_xml_node()])
        return result

    @classmethod
    def from_xml_node(cls, node):
        result = cls(node.get('id'), is_true(node.get('zeroCost')))
        n = get_unique_child(node, MachineryComposition.tag, False)
        if n is not None:
            result.machinery_composition = MachineryComposition.from_xml_node(n)
        n = get_unique_child(node, EnzymaticActivity.tag)
        result.enzymatic_activity = EnzymaticActivity.from_xml_node(n)
        return result

class ListOfEnzymes(ListOf):
    tag = 'listOfEnzymes'
    list_element = Enzyme

class EnzymaticActivity(object):
    tag = 'enzymaticActivity'
    
    def __init__(self, reaction):
        self.reaction = reaction
        self.enzyme_efficiencies = ListOfEnzymeEfficiencies()
        self.transporter_efficiency = TransporterEfficiency()

    def to_xml_node(self):
        result = etree.Element(self.tag)
        result.set('reaction', self.reaction)
        # optional subelements
        opt = [self.enzyme_efficiencies, self.transporter_efficiency]
        result.extend([e.to_xml_node() for e in opt if not(e.is_empty())])
        return result

    @classmethod
    def from_xml_node(cls, node):
        reaction = node.get('reaction')
        result = cls(reaction)
        n = get_unique_child(node, ListOfEnzymeEfficiencies.tag)
        result.enzyme_efficiencies = ListOfEnzymeEfficiencies.from_xml_node(n)
        n = get_unique_child(node, TransporterEfficiency.tag, False)
        if n is not None: result.transporter_efficiency = \
           TransporterEfficiency.from_xml_node(n)
        return result

class EnzymeEfficiency(object):
    tag = 'enzymeEfficiency'
    
    def __init__(self, function, parameters = {}):
        self.function = function
        self.parameters = ListOfParameters()
        for key,value in parameters.items():
            self.parameters.append(Parameter(key, value))

    def to_xml_node(self):
        result = etree.Element(self.tag)
        result.set('function', self.function)
        result.extend([self.parameters.to_xml_node()])
        return result

    @classmethod
    def from_xml_node(cls, node):
        result = cls(node.get('function'))
        n = get_unique_child(node, ListOfParameters.tag)
        result.parameters = ListOfParameters.from_xml_node(n)
        return result
        
class ListOfEnzymeEfficiencies(ListOf):
    tag = 'listOfEnzymeEfficiencies'
    list_element = EnzymeEfficiency

class TransporterEfficiency(ListOf):
    tag = 'transporterEfficiency'
    list_element = Function
