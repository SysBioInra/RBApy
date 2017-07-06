"""
Module containing valid RBA functions and Functions container class.
"""

# python 2/3 compatibility
from __future__ import division, print_function

# global imports
import numpy

class Functions(object):
    """
    Class building and storing RBA functions (including aggregates).

    Attributes:
        functions: dict mapping function or aggregate id with object used to
            compute it.
    """
    def __init__(self, functions, aggregates):
        """
        Constructor.

        Args:
            functions: xml structure containing function information.
            aggregates: xml structure containing aggregate information.
        """
        self.functions = {}
        for fn in functions:
            params = {p.id: p.value for p in fn.parameters}
            self.functions[fn.id] = build_function(fn.type, params, fn.variable)
        for fn in aggregates:
            self.functions[fn.id] = build_aggregate(fn, self.functions)

    def __getitem__(self, fn_id):
        """
        Get function or aggregate matching given id.

        Args:
            fn_id: id of function to retrieve.

        Returns:
            Function object.
        """
        return self.functions[fn_id]

class ConstantFunction(object):
    """
    Class computing constant functions.

    Attributes:
        name: identifier of this class.
        variable: variable used by function (if applicable).
    """

    name = 'constant'

    def __init__(self, parameters, variable=None):
        """
        Constructor.

        Args:
            parameters: dict used to initialize function. It must contain the
                the following keys: CONSTANT.
            variable: variable used by function (if applicable).
        """
        self.variable = variable
        self._constant = parameters['CONSTANT']

    def evaluate(self, x):
        """
        Evaluate function at given point.

        Args:
            x: point where function is to be computed.

        Returns:
            Value of the function at given point.
        """
        return self._constant

class ExponentialFunction(object):
    """
    Class computing exponential functions.

    Attributes:
        name: identifier of this class.
        variable: variable used by function (if applicable).
    """

    name = 'exponential'

    def __init__(self, parameters, variable=None):
        """
        Constructor.

        Args:
            parameters: dict used to initialize function. It must contain the
                the following keys: RATE.
            variable: variable used by function (if applicable).
        """
        self.variable = variable
        self._rate = parameters['RATE']

    def evaluate(self, x):
        """
        Evaluate function at given point.

        Args:
            x: point where function is to be computed.

        Returns:
            Value of the function at given point.
        """
        return numpy.exp(self._rate*x)

class IndicatorFunction(object):
    """
    Class computing indicator functions.

    Attributes:
        name: identifier of this class.
        variable: variable used by function (if applicable).
    """

    name = 'indicator'

    def __init__(self, parameters, variable=None):
        """
        Constructor.

        Args:
            parameters: dict used to initialize function. It must contain the
                the following keys: X_MIN, X_MAX.
            variable: variable used by function (if applicable).
        """
        self.variable = variable
        self._x_min = parameters['X_MIN']
        self._x_max = parameters['X_MAX']

    def evaluate(self, x):
        """
        Evaluate function at given point.

        Args:
            x: point where function is to be computed.

        Returns:
            Value of the function at given point.
        """
        return (x > self._x_min) and (x < self._x_max)

class LinearFunction(object):
    """
    Class computing linear functions.

    Attributes:
        name: identifier of this class.
        variable: variable used by function (if applicable).
    """

    name = 'linear'

    def __init__(self, parameters, variable=None):
        """
        Constructor.

        Args:
            parameters: dict used to initialize function. It must contain the
                the following keys: X_MIN, X_MAX, LINEAR_COEF, LINEAR_CONSTANT,
                Y_MIN, Y_MAX.
            variable: variable used by function (if applicable).
        """
        self.variable = variable
        self._x_min = parameters['X_MIN']
        self._x_max = parameters['X_MAX']
        self._coef = parameters['LINEAR_COEF']
        self._constant = parameters['LINEAR_CONSTANT']
        self._y_min = parameters['Y_MIN']
        self._y_max = parameters['Y_MAX']

    def evaluate(self, x):
        """
        Evaluate function at given point.

        Args:
            x: point where function is to be computed.

        Returns:
            Value of the function at given point.
        """
        x_eval = min(max(x, self._x_min), self._x_max)
        y = self._coef * x_eval + self._constant
        return min(max(y, self._y_min), self._y_max)

class MichaelisMentenFunction(object):
    """
    Class computing michaelis menten functions.

    Attributes:
        name: identifier of this class.
        variable: variable used by function (if applicable).
    """

    name = 'michaelisMenten'

    def __init__(self, parameters, variable=None):
        """
        Constructor.

        Args:
            parameters: dict used to initialize function. It must contain the
                the following keys: kmax, Km. Optionally, it may contain Y_MIN.
            variable: variable used by function (if applicable).
        """
        self.variable = variable
        self._kmax = float(parameters['kmax'])
        self._km = float(parameters['Km'])
        self._y_min = parameters.get('Y_MIN', None)

    def evaluate(self, x):
        """
        Evaluate function at given point.

        Args:
            x: point where function is to be computed.

        Returns:
            Value of the function at given point.
        """
        y = self._kmax * x / (x + self._km)
        return max(y, self._y_min) if self._y_min else y

class MultiplicationFunction(object):
    """
    Class computing multiplication functions.

    Attributes:
        name: identifier of this class.
        variable: variable used by function (if applicable).
    """

    name = 'multiplication'

    def __init__(self, function_handles):
        """
        Constructor.

        Args:
            function_handles: list of functions that have to be multiplied.
        """
        self.variable = None
        self._operands = function_handles

    def evaluate(self, x):
        """
        Evaluate function at given point.

        Args:
            x: point where function is to be computed.

        Returns:
            Value of the function at given point.
        """
        y = 1
        for op in self._operands:
            y *= op.evaluate(x)
        return y

# list of accepted function names and classes implementing them
SIMPLE = [ConstantFunction, LinearFunction, IndicatorFunction,
          ExponentialFunction, MichaelisMentenFunction]
AGGREGATE = [MultiplicationFunction]
VALID_FNS = {c.name: c for c in SIMPLE}
VALID_AGGS = {c.name: c for c in AGGREGATE}

def build_function(type_, params, variable=None):
    """
    Create object function matching type and parameters given.

    Args:
        type_: string matching name attribute of one of the classes of this
             module.
        params: dict mapping parameter names with their values.
        variable: name of variable used by function (if applicable).

    Returns:
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
        print('Missing parameter: ' + error.message)
        raise UserWarning('Invalid function.')

def build_aggregate(agg, known_functions):
    """
    Create aggregate from xml_structure.

    Args:
        agg: xml structure describing aggregate.
        known_functions: dict mapping known function ids with object used to
            compute them. These functions are used to build the aggregate.

    Returns:
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
        fn_handles = [known_functions[ref.function] \
                      for ref in agg.function_references]
    except KeyError as error:
        print('Unknown function: ' + error.args[0])
        raise UserWarning('Invalid aggregate.')
    return agg_class(fn_handles)
