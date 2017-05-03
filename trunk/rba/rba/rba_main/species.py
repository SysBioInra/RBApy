
from scipy.sparse import csr_matrix, lil_matrix, coo_matrix, hstack, eye
import numpy

class Machinery(object):
    def __init__(self, composition, processing_cost, weight):
        self.composition = composition
        self.processing_cost = processing_cost
        self.weight = weight

    def add_null_machineries(self, number):
        self.composition = self._pad_matrix(self.composition, number)
        self.processing_cost = self._pad_matrix(self.processing_cost, number)
        self.weight = self._pad_matrix(self.weight, number)

    @staticmethod
    def _pad_matrix(M, nb_cols):
        return hstack([M, coo_matrix((M.shape[0], nb_cols))]).tocsr()


class Species(object):
    def __init__(self, data, metabolites):
        self._metabolites = metabolites
        self._compartments = [c.id for c in data.metabolism.compartments]
        # extract composition of base species
        self.ids = metabolites \
                   + [m.id for m in data.proteins.macromolecules] \
                   + [m.id for m in data.rnas.macromolecules] \
                   + [m.id for m in data.dna.macromolecules]
        # metabolites
        species_comp = -eye(len(metabolites))
        species_proc = csr_matrix((len(data.processes.processes),
                                   len(metabolites)))
        species_deg = eye(len(metabolites))
        species_deg_proc = csr_matrix((len(data.processes.processes),
                                       len(metabolites)))
        species_weight = csr_matrix((len(self._compartments),
                                     len(metabolites)))
        # macromolecules
        [macro_comp, macro_proc, macro_deg, macro_deg_proc, macro_weight] \
            = self._macromolecule_composition(data)
        self.production = hstack([species_comp, macro_comp]).tocsr()
        self.prod_proc_cost = hstack([species_proc, macro_proc]).tocsr()
        self.degradation = hstack([species_deg, macro_deg]).tocsr()
        self.deg_proc_cost = hstack([species_deg_proc, macro_deg_proc]).tocsr()
        self.weight = hstack([species_weight, macro_weight]).tocsr()

    def create_machinery(self, machinery_set):
        species = lil_matrix((len(self.ids), len(machinery_set)))
        for col, m in enumerate(machinery_set):
            for r in m.reactants:
                species[self.ids.index(r.species), col] += r.stoichiometry
            for p in m.products:
                species[self.ids.index(p.species), col] -= p.stoichiometry
        return Machinery(self.production*species,
                         self.prod_proc_cost*species,
                         self.weight*species)

    def metabolite_synthesis(self):
        """
        """
        names = []
        reactions = []
        for m in self._metabolites:
            # if a metabolite is also a macromolecule, its id will appear
            # twice in the species list
            indices = [ind for (ind,i) in enumerate(self.ids) if i == m]
            if len(indices) == 1: continue
            # create biosynthesis reaction
            reaction = self.production[:, indices[1]].tolil()
            reaction[indices[0],0] = 1
            reactions.append(reaction)
            names.append(m + '_synthesis')
        return (reactions, names)

    def _component_maps(self, components, data):
        component_map = {}
        for m in data.processes.component_maps: component_map[m.id] = m
        production_map = {}
        degradation_map = {}
        nb_processes = len(data.processes.processes)
        for s in components:
            production_map[s] = ComponentMap(components[s], self._metabolites,
                                             nb_processes)
            degradation_map[s] = ComponentMap(components[s], self._metabolites,
                                              nb_processes)
        for p_index, process in enumerate(data.processes.processes):
            for n in process.operations.productions:
                production_map[n.set].add(component_map[n.component_map],
                                          p_index)
            for n in process.operations.degradations:
                degradation_map[n.set].add(component_map[n.component_map],
                                           p_index)            
        return (production_map, degradation_map)
        
    def _macromolecule_components(self, macro_set, data):
        # useful information
        nb_processes = len(data.processes.processes)
        components = [c.id for c in macro_set.components]
        # extract component matrix and location information
        nb_macros = len(macro_set.macromolecules)
        C = lil_matrix((len(components), nb_macros))
        compartment = []
        for col, m in enumerate(macro_set.macromolecules):
            for c in m.composition:
                C[components.index(c.component), col] = c.stoichiometry
            compartment.append(self._compartments.index(m.compartment))
        C = C.tocsr()
        # compute weight and associate weight with location
        w = csr_matrix([c.weight for c in macro_set.components], dtype='float')
        W = csr_matrix(((w*C).toarray().ravel(),
                        (compartment, range(nb_macros))),
                       shape=(len(self._compartments), nb_macros))
        return C, W

    def _macromolecule_composition(self, data):
        # get base matrices
        components = {}
        component_names = {}
        weights = {}
        subdata = {'protein': data.proteins, 'rna': data.rnas, 'dna': data.dna}
        sets = ['protein', 'rna', 'dna']
        for s in sets:
            [components[s], weights[s]] \
                = self._macromolecule_components(subdata[s], data)
            component_names[s] = [c.id for c in subdata[s].components]
        [production_map, degradation_map] \
            = self._component_maps(component_names, data)
        # compute composition and weight
        production_metabolites = []
        production_cost = []
        degradation_metabolites = []
        degradation_cost = []
        weight = []
        for s in sets:
            [comp, cost] = production_map[s].apply_map(components[s])
            production_metabolites.append(comp)
            production_cost.append(cost)
            [comp, cost] = degradation_map[s].apply_map(components[s])
            degradation_metabolites.append(comp)
            degradation_cost.append(cost)
            weight.append(weights[s])
        return (hstack(production_metabolites), hstack(production_cost),
                hstack(degradation_metabolites), hstack(degradation_cost),
                hstack(weight))

class ComponentMap(object):
    """
    Class used to store component maps.
    """
    def __init__(self, components, metabolites, nb_processes):
        """
        Constructor from XML node.
        :param xml_node: ComponentMap XML node in RBA format.
        :param metabolites: list of existing metabolites.
        :type xml_node: xml.dom.minidom element.
        :type metabolites: list of string
        """
        nb_metabolites = len(metabolites)
        nb_components = len(components)
        self._metabolites = metabolites
        self._metabolite_constant = numpy.zeros(nb_metabolites)
        self._processing_constant = numpy.zeros(nb_processes)
        self._components = components
        self._metabolite_table = numpy.zeros([nb_metabolites, nb_components])
        self._processing_table = numpy.zeros([nb_processes, nb_components])

    def add(self, map_, process_index):
        # store constant costs
        self._metabolite_constant += self._cost_vector(map_.constant_cost)
        # store component based costs
        for c in map_.costs:
            c_index = self._components.index(c.component)
            self._processing_table[process_index, c_index] += c.processing_cost
            self._metabolite_table[:, c_index] += self._cost_vector(c)

    def _cost_vector(self, cost):
        result = numpy.zeros(len(self._metabolites))
        for r in cost.reactants:
            result[self._metabolites.index(r.species)] -= r.stoichiometry
        for r in cost.products:
            result[self._metabolites.index(r.species)] += r.stoichiometry
        return result

    def apply_map(self, component_matrix):
        """
        """
        col_selector = numpy.zeros(component_matrix.shape[1])
        metab_cost = csr_matrix(self._metabolite_table) * component_matrix \
                     + csr_matrix(self._metabolite_constant).T[:, col_selector]
        proc_cost = csr_matrix(self._processing_table) * component_matrix \
                    + csr_matrix(self._processing_constant).T[:, col_selector]
        return [metab_cost, proc_cost]
