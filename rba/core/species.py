"""
Module defining Machinery, Species and ComponentMap classes.
"""

# python 2/3 compatibility
from __future__ import division, print_function

# global imports
from collections import namedtuple
from scipy.sparse import csr_matrix, lil_matrix, hstack, eye
import numpy

# class storing machinery-related information
Machinery = namedtuple('Machinery', 'composition processing_cost weight')

class Species(object):
    """
    Class storing species-related information.

    Attributes:
        ids: identifiers of species stored (metabolites/macromolecules).
        production: production matrix (in terms of metabolites).
        prod_proc_cost: production processing cost matrix.
        degradation: degradation matrix (in terms of metabolites).
        deg_proc_cost: degradation processing cost matrix.
        weight: weight matrix.
    """

    def __init__(self, data, metabolites):
        self._metabolites = metabolites
        self._compartments = [c.id for c in data.metabolism.compartments]
        # extract composition of base species
        self.ids = (metabolites
                    + [m.id for m in data.proteins.macromolecules]
                    + [m.id for m in data.rnas.macromolecules]
                    + [m.id for m in data.dna.macromolecules])
        # polymers and metabolites are allowed to have the same identifier
        # by looping on reversed list, we ensure that the index of the
        # metabolite is returned
        self._index = {m: i for i, m in reversed(list(enumerate(self.ids)))}
        # metabolites (weights and processing costs are zero)
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
        """
        Create machinery from list of RBA machinery composition structures.

        Args:
            machinery_set: list of RBA machinery composition structures.

        Returns:
            Machinery object containing composition, processing cost and
            weight matrices.
        """
        species = lil_matrix((len(self.ids), len(machinery_set)))
        for col, machinery in enumerate(machinery_set):
            for reac in machinery.reactants:
                species[self._index[reac.species], col] += reac.stoichiometry
            for prod in machinery.products:
                species[self._index[prod.species], col] -= prod.stoichiometry
        return Machinery(self.production*species,
                         self.prod_proc_cost*species,
                         self.weight*species)

    def metabolite_synthesis(self):
        """
        Create reactions corresponding to synthesis of macrometabolites.

        Macrometabolites are species that are both a metabolite and a
        macromolecule (typically tRNAs).

        Returns:
            Tuple where first element is a list of stoichiometries vectors,
            each vector representing a reaction. The second element are the
            ids of the metabolites being synthesized by these reactions.
        """
        names = []
        reactions = []
        nb_met = len(self._metabolites)
        macrometabolites = self.ids[nb_met:]
        for index, macro in enumerate(macrometabolites):
            # if a macromolecule is also a metabolite,
            # it appears twice in the species list.
            met_index = self._index[macro]
            macro_index = nb_met + index
            if met_index < nb_met:
                # create biosynthesis reaction
                reaction = self.production[:, macro_index].tolil()
                reaction[met_index, 0] = 1
                reactions.append(reaction)
                names.append(macro + '_synthesis')
        return (reactions, names)

    def _component_maps(self, component_sets, processes):
        """
        Convert component map data to ComponentMap objects.
        """
        component_map = {m.id: m for m in processes.component_maps}
        production_map = {}
        degradation_map = {}
        nb_processes = len(processes.processes)
        for set_, components in component_sets.items():
            production_map[set_] = ComponentMap(components, self._metabolites,
                                                nb_processes)
            degradation_map[set_] = ComponentMap(components, self._metabolites,
                                                 nb_processes)
        for p_index, process in enumerate(processes.processes):
            for prod in process.operations.productions:
                production_map[prod.set].add(component_map[prod.component_map],
                                             p_index)
            for deg in process.operations.degradations:
                degradation_map[deg.set].add(component_map[deg.component_map],
                                             p_index)
        return (production_map, degradation_map)

    def _macromolecule_components(self, macro_set):
        """
        Extract component and weight matrices from macromolecule data.

        A component matrix is the description of macromolecules in terms
        of components (e.g. amino acids). Compare composition matrix, the
        description  in terms of metabolites consumed and
        produced for synthesizing one macromolecule.
        """
        # useful information
        components = [c.id for c in macro_set.components]
        # extract component matrix and location information
        nb_macros = len(macro_set.macromolecules)
        C = lil_matrix((len(components), nb_macros))
        compartment = []
        for col, macro in enumerate(macro_set.macromolecules):
            for c in macro.composition:
                C[components.index(c.component), col] = c.stoichiometry
            compartment.append(self._compartments.index(macro.compartment))
        C = C.tocsr()
        # compute weight and associate weight with location
        # we first compute weight per component, then weight per molecule
        w = csr_matrix([c.weight for c in macro_set.components], dtype='float')
        W = csr_matrix(((w*C).toarray().ravel(),
                        (compartment, range(nb_macros))),
                       shape=(len(self._compartments), nb_macros))
        return C, W

    def _macromolecule_composition(self, data):
        """
        Compute production, degradation, processing_costs and weight of
        macromolecules.
        """
        # get base matrices
        components = {}
        component_names = {}
        weights = {}
        subdata = {'protein': data.proteins, 'rna': data.rnas, 'dna': data.dna}
        sets = ['protein', 'rna', 'dna']
        for set_ in sets:
            [components[set_], weights[set_]] \
                = self._macromolecule_components(subdata[set_])
            component_names[set_] = [c.id for c in subdata[set_].components]
        [production_map, degradation_map] \
            = self._component_maps(component_names, data.processes)
        # compute composition and weight
        production_metabolites = []
        production_cost = []
        degradation_metabolites = []
        degradation_cost = []
        weight = []
        for set_ in sets:
            [comp, cost] = production_map[set_].apply_map(components[set_])
            production_metabolites.append(comp)
            production_cost.append(cost)
            [comp, cost] = degradation_map[set_].apply_map(components[set_])
            degradation_metabolites.append(comp)
            degradation_cost.append(cost)
            weight.append(weights[set_])
        return (hstack(production_metabolites), hstack(production_cost),
                hstack(degradation_metabolites), hstack(degradation_cost),
                hstack(weight))


class ComponentMap(object):
    """
    Class storing component maps.
    """

    def __init__(self, components, metabolites, nb_processes):
        """
        Constructor.

        Args:
            components: list of components handled by component map.
            metabolites: list of existing metabolites.
            nb_processes: number of processes.
        """
        nb_metabolites = len(metabolites)
        nb_components = len(components)
        self._met_index = {m: i for i, m in enumerate(metabolites)}
        self._metabolite_constant = numpy.zeros(nb_metabolites)
        self._processing_constant = numpy.zeros(nb_processes)
        self._components = components
        self._metabolite_table = numpy.zeros([nb_metabolites, nb_components])
        self._processing_table = numpy.zeros([nb_processes, nb_components])

    def add(self, map_, process_index):
        """
        Add item to component map.

        Args:
            map_: component map as stored in RBA data.
            process_index: index of process using this component map.
        """
        # store constant costs
        self._metabolite_constant += self._cost_vector(map_.constant_cost)
        # store component based costs
        for cost in map_.costs:
            c_index = self._components.index(cost.component)
            self._processing_table[process_index, c_index] \
                += cost.processing_cost
            self._metabolite_table[:, c_index] += self._cost_vector(cost)

    def _cost_vector(self, cost):
        """
        Transform cost data into a metabolite vector.
        """
        result = numpy.zeros(len(self._met_index))
        for reac in cost.reactants:
            result[self._met_index[reac.species]] -= reac.stoichiometry
        for prod in cost.products:
            result[self._met_index[prod.species]] += prod.stoichiometry
        return result

    def apply_map(self, component_matrix):
        """
        Transform component matrix to metabolite matrix.

        Args:
            component_matrix: description of macromolecules in terms of
                components (columns are macromolecules, rows are components).

        Returns:
            Tuple (composition, processing_cost) where:
            composition is a metabolite matrix
            describing metabolites consumed/produced during macromolecule
            synthesis/degradation (depending on definition of the map).
            Columns are macromolecules. Rows are metabolites. A negative
            coefficient means the metabolite is *produced* (this is a
            composition matrix, not a reaction matrix).
            processing_cost is a matrix where columns are macromolecules and
            lines are processes. It describes how many resources of a process
            are used during macromolecule synthesis/degradation.
        """
        # column selector used to duplicate vectors to match final matrix size
        cols = numpy.zeros(component_matrix.shape[1])
        metab_cost = (csr_matrix(self._metabolite_table) * component_matrix
                      + csr_matrix(self._metabolite_constant).T[:, cols])
        proc_cost = (csr_matrix(self._processing_table) * component_matrix
                     + csr_matrix(self._processing_constant).T[:, cols])
        return [metab_cost, proc_cost]
