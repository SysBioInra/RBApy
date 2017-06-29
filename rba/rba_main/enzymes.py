"""
Module defining Enzymes and EnzymeEfficiency classes.
"""

# python 2/3 compatibility
from __future__ import division, print_function

# global imports
import numpy

# local imports
from rba.rba_main.functions import build_function

class Enzymes(object):
    """
    Class computing enzyme-related substructures.

    Attributes:
        ids: list of identifiers of enzymes having a nonzero production cost.
        reaction_catalyzed: list of identifiers of reaction catalyzed by enzymes
            in ids.
        machinery (species.Machinery): object containing composition-related
            information of enzymes.
        efficiency (EnzymeEfficiency): object used to compute enzyme
            efficiencies depending on medium and growth-rate.
    """
    def __init__(self, enzymes, species, reactions):
        """
        Constructor.

        Args:
            enzymes: xml structure containing enzyme information.
            species (species.Species): object containing all information about
                chemical species in the system.
            reactions: list of reaction identifiers.
        """
        # check that all reactions are found and keep only enzymes
        # that have a machinery
        reactions_left = reactions[:]
        nonzero_enzymes = []
        self.reaction_catalyzed = []
        for enzyme in enzymes.enzymes:
            reaction = enzyme.enzymatic_activity.reaction
            reactions_left.remove(reaction)
            if not(enzyme.zero_cost or enzyme.machinery_composition.is_empty()):
                nonzero_enzymes.append(enzyme)
                self.reaction_catalyzed.append(reaction)
        self.ids = [e.id for e in nonzero_enzymes]
        if reactions_left:
            print('Warning: did not find enzymes for following reactions: '
                  + ', '.join(reactions_left))

        # extract machinery information
        machinery = [e.machinery_composition for e in nonzero_enzymes]
        self.machinery = species.create_machinery(machinery)

        # extract efficiency information
        self.efficiency = EnzymeEfficiency(enzymes.efficiency_functions)
        for enzyme in nonzero_enzymes:
            self.efficiency.add_activity(enzyme.enzymatic_activity)

class EnzymeEfficiency(object):
    """
    Class computing efficiency of enzymes depending on medium and growth rate.
    """

    def __init__(self, efficiency_functions):
        """
        Constructor.

        Args:
            efficiency_functions: list of functions.Function describing ids and
                types of all efficiency functions that an enzyme must define.
        """
        # read list of efficiency functions
        self._fn_types = {fn.id: fn.type for fn in efficiency_functions}
        self._efficiency = {fn_id: [] for fn_id in self._fn_types}
        self._import = []
        # use default values for import and efficiency function
        self._eff_fn = self._efficiency[efficiency_functions[0].id]
        self._import_values = 1

    def add_activity(self, enzymatic_activity):
        """
        Add enzymatic activity.

        Args:
            enzymatic_activity: xml structure containing enzymatic activity of
                an enzyme. The structure must define parameters for all
                efficiency functions provided at construction.
        """
        # base efficiency
        params = {}
        for fn in enzymatic_activity.enzyme_efficiencies:
            params[fn.function] = {p.id: p.value for p in fn.parameters}
        for fn_id in self._efficiency:
            try:
                self._efficiency[fn_id].append \
                    (build_function(self._fn_types[fn_id], params[fn_id]))
            except KeyError:
                print('Missing parameters for enzymatic function ' + fn_id)
                raise UserWarning('Invalid enzyme file.')
        # import efficiency
        t_eff = enzymatic_activity.transporter_efficiency
        if t_eff:
            fns = []
            for fn in t_eff:
                params = {p.id: p.value for p in fn.parameters}
                fns.append(build_function(fn.type, params, fn.variable))
            self._import.append(fns)
        else:
            self._import.append([])

    def set_function(self, fn_id):
        """
        Set efficiency function that should be used for computations.

        Args:
            fn_id: id of efficiency function that must match one of the
                functions defined at construction.
        """
        self._eff_fn = self._efficiency[fn_id]

    def update_import(self, concentration):
        """
        Update transport terms based on medium concentrations.

        Args:
            concentration: dict mapping metabolite prefixes with their
                concentration.
        """
        self._import_values = numpy.ones(len(self._import))
        for i, import_fn in enumerate(self._import):
            for i_fn in import_fn:
                # /!\ we identify metabolites by their prefix !!!
                key = i_fn.variable.rsplit('_', 1)[0]
                self._import_values[i] *= i_fn.evaluate(concentration[key])

    def compute(self, mu):
        """
        Compute efficiency for given growth rate.

        Args:
            mu: growth rate

        Returns:
            Vector containing enzyme efficiencies in the order they have been
            added.
        """
        efficiency = numpy.ones(len(self._eff_fn))
        for i, fn in enumerate(self._eff_fn):
            if fn:
                efficiency[i] = fn.evaluate(mu)
        return (efficiency*self._import_values, efficiency)
