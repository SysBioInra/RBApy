"""Module with parameter container class."""

# python 2/3 compatibility
from __future__ import division, print_function, absolute_import

# local imports
from rba.core.functions import build_function, build_aggregate


class Parameters(object):
    """
    Class storing RBA parameters (including aggregates).

    Attributes
    ----------
    parameters : dict
        Mapping from parameter id with object used to compute it.

    """

    def __init__(self, functions, aggregates):
        """
        Constructor.

        Parameters
        ----------
        functions : rba.xml.ListOfFunctions
            structure containing function information.
        aggregates: rba.xml.ListOfAggregates
            structure containing aggregate information.

        """
        self.parameters = {}
        self._growth_rate_fn = []
        self._medium_fn = []
        self._growth_rate_agg = []
        self._medium_agg = []
        for fn in functions:
            params = {p.id: p.value for p in fn.parameters}
            new_fn = build_function(fn.type, params, fn.variable)
            self.parameters[fn.id] = new_fn
            if new_fn.is_growth_rate_dependent():
                self._growth_rate_fn.append(new_fn)
            if new_fn.is_medium_dependent():
                self._medium_fn.append(new_fn)
        for fn in aggregates:
            new_agg = build_aggregate(fn, self.parameters)
            self.parameters[fn.id] = new_agg
            if new_agg.is_growth_rate_dependent():
                self._growth_rate_agg.append(new_agg)
            if new_agg.is_medium_dependent():
                self._medium_agg.append(new_agg)

    def __getitem__(self, parameter_id):
        """
        Get function or aggregate matching given id.

        Parameters
        ----------
        fn_id : str
            id of parameter to retrieve.

        Returns
        -------
        Function object.

        """
        return self.parameters[parameter_id]

    def update_growth_rate(self, growth_rate):
        """
        Compute parameters for given growth rate.

        Parameters
        ----------
        growth_rate : float
            Current growth rate.

        """
        # update functions first, then aggregates !!!
        for fn in self._growth_rate_fn:
            fn.update(growth_rate)
        for agg in self._growth_rate_agg:
            agg.update()

    def update_medium(self, medium):
        """
        Compute parameters for given medium concentrations.

        Parameters
        ----------
        concentration : dict
            Mapping of metabolite prefixes with their concentration.

        """
        # update functions first, then aggregates !!!
        for fn in self._medium_fn:
            # /!\ we identify metabolites by their prefix !!!
            prefix = fn.variable.rsplit('_', 1)[0]
            fn.update(medium[prefix])
        for agg in self._medium_agg:
            agg.update()
