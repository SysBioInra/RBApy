"""Module defining Machinery, Species and ProcessingMap classes."""

# python 2/3 compatibility
from __future__ import division, print_function, absolute_import

# global imports
from collections import namedtuple
from scipy.sparse import (
    csr_matrix, lil_matrix, coo_matrix, hstack, eye
    )
import numpy

# class storing machinery-related information
Machinery = namedtuple('Machinery', 'composition processing_cost weight')


class Species(object):
    """
    Species-related information.

    Parameters
    ----------
    ids : list of str
        Identifiers of species stored (metabolites/macromolecules).
    production : sparse matrix
        Production matrix (in terms of metabolites).
    prod_proc_cost : sparse matrix
        Production processing cost matrix.
    degradation : sparse matrix
        Degradation matrix (in terms of metabolites).
    deg_proc_cost : sparsematrix
        Degradation processing cost matrix.
    weight : sparse matrix
        Weight matrix (in terms of compartments).

    """

    def __init__(self, data, metabolites):
        """Constructor."""
        self._metabolites = metabolites
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
        nb_comp = len(data.metabolism.compartments)
        nb_met = len(metabolites)
        nb_processes = len(data.processes.processes)
        met_comp = -eye(nb_met)
        met_proc = csr_matrix((nb_processes, nb_met))
        met_deg = eye(nb_met)
        met_deg_proc = csr_matrix((nb_processes, nb_met))
        met_weight = csr_matrix((nb_comp, nb_met))
        # macromolecules
        [macro_comp, macro_proc, macro_deg, macro_deg_proc, macro_weight] \
            = compute_macromolecule_composition(data, metabolites)
        self.production = hstack([met_comp, macro_comp]).tocsr()
        self.prod_proc_cost = hstack([met_proc, macro_proc]).tocsr()
        self.degradation = hstack([met_deg, macro_deg]).tocsr()
        self.deg_proc_cost = hstack([met_deg_proc, macro_deg_proc]).tocsr()
        self.weight = hstack([met_weight, macro_weight]).tocsr()

    def create_machinery(self, machinery_set):
        """
        Create machineries from a list of RBA machinery composition structures.

        Parameters
        ----------
        machinery_set : list of rba.xml.MachineryComposition
            Machinery compositions.

        Returns
        -------
        Machinery object
            Contains the composition, processing cost and weight matrices of
            machineries provided as input.

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

        Returns
        -------
        Tuple of 2 elements
            First element is a list of stoichiometries vectors,
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
        return reactions, names


def compute_macromolecule_composition(data, metabolites):
    """
    Compute base information of macromolecules.

    Returns
    -------
    (production, production_processing_cost, degradation,
     degradation_processing_cost, weight) tuple

    """
    nb_processes = len(data.processes.processes)
    compartments = [c.id for c in data.metabolism.compartments]
    # get base macromolecule information
    proteins = MacromoleculeSet(data.proteins, compartments,
                                metabolites, nb_processes)
    rnas = MacromoleculeSet(data.rnas, compartments, metabolites, nb_processes)
    dna = MacromoleculeSet(data.dna, compartments, metabolites, nb_processes)
    # apply processing maps
    macro_sets = {'protein': proteins, 'rna': rnas, 'dna': dna}
    maps = {m.id: m for m in data.processes.processing_maps}
    for p_index, process in enumerate(data.processes.processes):
        for prod in process.processings.productions:
            inputs = [i.species for i in prod.inputs]
            macro_sets[prod.set].apply_production_map(
                maps[prod.processing_map], p_index, inputs
                )
        for deg in process.processings.degradations:
            inputs = [i.species for i in deg.inputs]
            macro_sets[deg.set].apply_degradation_map(
                maps[deg.processing_map], p_index, inputs
                )
    # aggregate matrices across sets
    production_metabolites = [s.production for s in (proteins, rnas, dna)]
    production_cost = [s.production_cost for s in (proteins, rnas, dna)]
    degradation_metabolites = [s.degradation for s in (proteins, rnas, dna)]
    degradation_cost = [s.degradation_cost for s in (proteins, rnas, dna)]
    weight = [s.weight for s in (proteins, rnas, dna)]
    return (hstack(production_metabolites), hstack(production_cost),
            hstack(degradation_metabolites), hstack(degradation_cost),
            hstack(weight))


class MacromoleculeSet(object):
    """Macromolecule information."""

    def __init__(self, macro_set, compartments, metabolites, nb_processes):
        """Initialize set with zero production/degradation costs."""
        self.components = [c.id for c in macro_set.components]
        self._molecule_index = {
            m.id: i for i, m in enumerate(macro_set.macromolecules)
            }
        self._component_matrix = self._extract_component_matrix(macro_set)
        self.weight = self._extract_weight_matrix(
            macro_set, self._component_matrix, compartments
            )
        self._metabolites = metabolites
        nb_met = len(metabolites)
        nb_mol = len(self._molecule_index)
        self.production = coo_matrix((nb_met, nb_mol))
        self.degradation = coo_matrix((nb_met, nb_mol))
        self.production_cost = coo_matrix((nb_processes, nb_mol))
        self.degradation_cost = coo_matrix((nb_processes, nb_mol))

    def apply_production_map(self, map_, process_index, inputs):
        self._apply_map(map_, inputs, process_index,
                        self.production, self.production_cost)

    def apply_degradation_map(self, map_, process_index, inputs):
        self._apply_map(map_, inputs, process_index,
                        self.degradation, self.degradation_cost)

    def _apply_map(self, map_, inputs, process_index, met_matrix, proc_matrix):
        # create column selector for inputs
        cols = numpy.array([self._molecule_index[i] for i in inputs],
                           dtype = int)
        proc_map = ProcessingMap(map_, self.components, self._metabolites)
        met, proc_cost = proc_map.apply_map(self._component_matrix[:, cols])
        # update production/degradation reactions
        met = met.tocoo()
        met_matrix.row = numpy.concatenate([met_matrix.row, met.row])
        met_matrix.col = numpy.concatenate([met_matrix.col, cols[met.col]])
        met_matrix.data = numpy.concatenate([met_matrix.data, met.data])
        # udpate procesing cost matrix
        if proc_cost.nnz:
            proc_cost = proc_cost.tocoo()
            proc_matrix.row = numpy.concatenate(
                [proc_matrix.row,
                 numpy.array([process_index]*len(proc_cost.data))]
                )
            proc_matrix.col = numpy.concatenate([proc_matrix.col,
                                                 cols[proc_cost.col]])
            proc_matrix.data = numpy.concatenate([proc_matrix.data,
                                                  proc_cost.data])

    def _extract_component_matrix(self, macro_set):
        """
        Extract component matrix from macromolecule data.

        A component matrix is the description of macromolecules in terms
        of components (e.g. amino acids). Compare composition matrix, the
        description  in terms of metabolites consumed and
        produced for synthesizing one macromolecule.

        """
        nb_macros = len(macro_set.macromolecules)
        C = lil_matrix((len(self.components), nb_macros))
        for col, macro in enumerate(macro_set.macromolecules):
            for c in macro.composition:
                C[self.components.index(c.component), col] = c.stoichiometry
        return C.tocsr()

    def _extract_weight_matrix(self, macro_set, C, compartments):
        """Compute weight and associate weight with location."""
        # we first compute weight per component, then weight per molecule
        w = csr_matrix([c.weight for c in macro_set.components], dtype='float')
        location = [compartments.index(m.compartment)
                    for m in macro_set.macromolecules]
        nb_macros = len(macro_set.macromolecules)
        W = csr_matrix(((w*C).toarray().ravel(),
                        (location, range(nb_macros))),
                       shape=(len(compartments), nb_macros))
        return W


class ProcessingMap(object):
    """Class storing processing maps."""

    def __init__(self, map_, components, metabolites):
        """
        Constructor.

        Parameters
        ----------
        map_ : rba.xml.ProcessingMap
            Structure containing processing map.
        components : list of rba.xml.Components
            Components handled by component map.
        metabolites : list of str
            Metabolites.

        """
        nb_metabolites = len(metabolites)
        nb_components = len(components)
        met_index = {m: i for i, m in enumerate(metabolites)}
        # store constant costs
        self._metabolite_constant \
            = self._cost_vector(map_.constant_processing, met_index)
        self._processing_constant = numpy.zeros(1)
        # store component based costs
        self._metabolite_table = numpy.zeros([nb_metabolites, nb_components])
        self._processing_table = numpy.zeros(nb_components)
        for proc in map_.component_processings:
            c_index = components.index(proc.component)
            self._processing_table[c_index] += proc.machinery_cost
            self._metabolite_table[:, c_index] += self._cost_vector(proc,
                                                                    met_index)

    def _cost_vector(self, proc, met_index):
        """Transform processing data into a metabolite vector."""
        result = numpy.zeros(len(met_index))
        for reac in proc.reactants:
            result[met_index[reac.species]] -= reac.stoichiometry
        for prod in proc.products:
            result[met_index[prod.species]] += prod.stoichiometry
        return result

    def apply_map(self, component_matrix):
        """
        Transform component matrix to metabolite matrix.

        Parameters
        ----------
        component_matrix: matrix
            Description of macromolecules in terms of components
            (columns are macromolecules, rows are components).

        Returns
        -------
        (composition, processing_cost) tuple
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
        cols = numpy.zeros(component_matrix.shape[1], dtype = int)
        metab_cost = (csr_matrix(self._metabolite_table) * component_matrix
                      + csr_matrix(self._metabolite_constant).T[:, cols])
        proc_cost = (csr_matrix(self._processing_table) * component_matrix
                     + csr_matrix(self._processing_constant).T[:, cols])
        return metab_cost, proc_cost
