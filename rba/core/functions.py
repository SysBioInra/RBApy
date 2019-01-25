"""Module containing valid RBA functions."""

# python 2/3 compatibility
from __future__ import division, print_function, absolute_import

# global imports
import numpy


class BaseFunction(object):
    """Mother of all functions."""

    def __init__(self, variable):
        """Build default object."""
        self.variable = variable

    def is_growth_rate_dependent(self):
        """Return whether function depends on growth rate."""
        return self.variable == 'growth_rate'

    def is_medium_dependent(self):
        """Return whether function depends on medium variable."""
        return self.variable and self.variable != 'growth_rate'


class ConstantFunction(BaseFunction):
    """
    Class computing constant functions.

    Attributes
    ----------
    name : str
        identifier of this class.
    variable : str
        variable used by function.
    value : float
        current function value.

    """

    name = 'constant'

    def __init__(self, parameters, variable=None):
        """
        Constructor.

        Parameters
        ----------
        parameters : dict
            Dict that must contain the following keys: CONSTANT.
        variable : str
            Function variable.

        """
        super(ConstantFunction, self).__init__(None)
        self.value = parameters['CONSTANT']

    def update(self, x):
        """Evaluate function."""
        pass


class ExponentialFunction(BaseFunction):
    """
    Class computing exponential functions.

    Attributes
    ----------
    name : str
        identifier of this class.
    variable : str
        variable used by function.
    value : float
        current function value.

    """

    name = 'exponential'

    def __init__(self, parameters, variable):
        """
        Constructor.

        Parameters
        ----------
        parameters : dict
            Dict that must contain the following keys: RATE.
        variable : str
            Function variable.

        """
        super(ExponentialFunction, self).__init__(variable)
        self.value = 1
        self._rate = parameters['RATE']

    def update(self, x):
        """Evaluate function."""
        self.value = numpy.exp(self._rate * x)


class IndicatorFunction(BaseFunction):
    """
    Class computing indicator functions.

    Attributes
    ----------
    name : str
        identifier of this class.
    variable : str
        variable used by function.
    value : bool
        current function value.

    """

    name = 'indicator'

    def __init__(self, parameters, variable):
        """
        Constructor.

        Parameters
        ----------
        parameters : dict
            Dict that must contain the following keys: X_MIN, X_MAX.
        variable : str
            Function variable.

        """
        super(IndicatorFunction, self).__init__(variable)
        self.value = 0
        self._x_min = parameters['X_MIN']
        self._x_max = parameters['X_MAX']

    def update(self, x):
        """Evaluate function."""
        self.value = (x > self._x_min) and (x < self._x_max)


class LinearFunction(BaseFunction):
    """
    Class computing linear functions.

    Attributes
    ----------
    name : str
        identifier of this class.
    variable : str
        variable used by function.
    value : float
        current function value.

    """

    name = 'linear'

    def __init__(self, parameters, variable):
        """
        Constructor.

        Parameters
        ----------
        parameters : dict
            Dict that must contain the following keys: X_MIN, X_MAX,
            LINEAR_COEF, LINEAR_CONSTANT, Y_MIN, Y_MAX.
        variable : str
            Function variable.

        """
        super(LinearFunction, self).__init__(variable)
        self.value = 0
        self._x_min = parameters['X_MIN']
        self._x_max = parameters['X_MAX']
        self._coef = parameters['LINEAR_COEF']
        self._constant = parameters['LINEAR_CONSTANT']
        self._y_min = parameters['Y_MIN']
        self._y_max = parameters['Y_MAX']

    def update(self, x):
        """Evaluate function."""
        x_eval = min(max(x, self._x_min), self._x_max)
        y = self._coef * x_eval + self._constant
        self.value = min(max(y, self._y_min), self._y_max)


class MichaelisMentenFunction(BaseFunction):
    """
    Class computing michaelis menten functions.

    Attributes
    ----------
    name : str
        identifier of this class.
    variable : str
        variable used by function.
    value : float
        current function value.

    """

    name = 'michaelisMenten'

    def __init__(self, parameters, variable):
        """
        Constructor.

        Parameters
        ----------
        parameters : dict
            Dict that must contain the following keys: kmax, Km.
            Optionally, it may contain Y_MIN.
        variable : str
            Function variable.

        """
        super(MichaelisMentenFunction, self).__init__(variable)
        self._kmax = float(parameters['kmax'])
        self._km = float(parameters['Km'])
        self._y_min = parameters.get('Y_MIN', None)
        self.value = 0

    def update(self, x):
        """Evaluate function at given point."""
        y = self._kmax * x / (x + self._km)
        self.value = max(y, self._y_min) if self._y_min else y


class CompetitiveInhibitionFunction(BaseFunction):
    """
    Class computing michaelis menten functions.

    Attributes
    ----------
    name : str
        identifier of this class.
    variable : str
        variable used by function.
    value : float
        current function value.

    """

    name = 'competitiveInhibition'

    def __init__(self, parameters, variable):
        """
        Constructor.

        Parameters
        ----------
        parameters : dict
            Dict that must contain the following keys: kmax, Km.
            Optionally, it may contain Y_MIN.
        variable : str
            Function variable.

        """
        super(CompetitiveInhibitionFunction, self).__init__(variable)
        self._kmax = float(parameters['kmax'])
        self._Km = float(parameters['Km'])
        self._Ki = float(parameters['Ki'])
        self._I = float(parameters['I'])
        self._y_min = parameters.get('Y_MIN', None)
        self.value = 0

    def update(self, x):
        """Evaluate function at given point."""
        y = self._kmax * x / (x + self._Km * (1 + self._I / self._Ki))
        self.value = max(y, self._y_min) if self._y_min else y


class InverseFunction(BaseFunction):
    """
    Class computing inverse functions.

    Attributes
    ----------
    name : str
        identifier of this class.
    variable : str
        variable used by function.
    value : float
        current function value.

    """

    name = 'inverse'

    def __init__(self, parameters, variable):
        """
        Constructor.

        Parameters
        ----------
        parameters : dict
            Dict that must contain the following keys: CONSTANT.
        variable : str
            Function variable.

        """
        super(InverseFunction, self).__init__(variable)
        self._xmax = parameters['CONSTANT']
        self.value = 0

    def update(self, x):
        """Evaluate function."""
        try:
            self.value = self._xmax / x
        except KeyError:
            print('variable is 0, impossible to do inversion')


class MultiplicationFunction(object):
    """
    Class computing multiplication functions.

    Attributes
    ----------
    name : str
        identifier of this class.
    value : int
        current function value.

    """

    name = 'multiplication'

    def __init__(self, function_handles):
        """
        Constructor.

        Parameters
        ----------
        function_handles : list of function handles
            Functions that have to be multiplied.

        """
        self._operands = function_handles
        self.value = 1
        self.update()

    def is_growth_rate_dependent(self):
        """Return whether function depends on growth rate."""
        return any(op.is_growth_rate_dependent() for op in self._operands)

    def is_medium_dependent(self):
        """Return whether function depends on medium variable."""
        return any(op.is_medium_dependent() for op in self._operands)

    def update(self):
        """Evaluate function."""
        y = 1
        for op in self._operands:
            y *= op.value
        self.value = y


# list of accepted function names and classes implementing them
SIMPLE = [ConstantFunction, LinearFunction, IndicatorFunction,
          ExponentialFunction, MichaelisMentenFunction, InverseFunction,
          CompetitiveInhibitionFunction]
AGGREGATE = [MultiplicationFunction]
VALID_FNS = {c.name: c for c in SIMPLE}
VALID_AGGS = {c.name: c for c in AGGREGATE}


def build_function(type_, params, variable):
    """
    Create object function matching type and parameters given.

    Parameters
    ----------
    type_ : str
        'name' attribute of one of the classes of this module.
    params : dict
        Dict mapping parameter names with their values.
    variable : str
        Function variable.

    Returns
    -------
    Function object matching type and parameters provided.

    """
    try:
        # retrieve class implementing function
        fn_class = VALID_FNS[type_]
    except KeyError:
        print('Unknown function type: ' + type_ + '. Valid types are: '
              + ', '.join(VALID_FNS.keys()))
        raise UserWarning('Invalid function.')
    try:
        return fn_class(params, variable)
    except KeyError as error:
        print('Missing parameter: ' + str(error))
        raise UserWarning('Invalid function.')


def build_aggregate(agg, known_functions):
    """
    Create aggregate from xml_structure.

    Parameters
    ----------
    agg : XML node
        Structure describing aggregate.
    known_functions : rba.core.parameters.Parameters
        Known parameters.

    Returns
    -------
    Aggregate object matching xml structure.

    """
    try:
        # retrieve class implementing aggregate
        agg_class = VALID_AGGS[agg.type]
    except KeyError:
        print('Unknown aggregate type: ' + agg.type + '. Valid types are: '
              + ', '.join(VALID_AGGS.keys()))
        raise UserWarning('Invalid aggregate.')
    try:
        # retrieve functions used in aggregate
        fn_handles = [known_functions[ref.function]
                      for ref in agg.function_references]
    except KeyError as error:
        print('Unknown function: ' + error.args[0])
        raise UserWarning('Invalid aggregate.')
    return agg_class(fn_handles)


zero_function = ConstantFunction({'CONSTANT': 0})
default_lb = ConstantFunction({'CONSTANT': -1e3})
default_ub = ConstantFunction({'CONSTANT': 1e5})
