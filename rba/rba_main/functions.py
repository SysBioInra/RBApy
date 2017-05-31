#
# File containing valid RBA functions.
#
import numpy
import abc

class Functions:
    def __init__(self, functions, aggregates):
        self.functions = {}
        for fn in functions:
            params = {p.id: p.value for p in fn.parameters}
            self.functions[fn.id] = build_function(fn.type, params, fn.variable)
        for fn in aggregates:
            self.functions[fn.id] = build_aggregate(fn, self.functions)
                                                          
    def __getitem__(self, fn_id):
        return self.functions[fn_id]

class RBAFunction(object):
    __metaclass__ = abc.ABCMeta

    def __init__(self, variable = None):
        self.variable = variable

    @abc.abstractmethod
    def evaluate(self, x): pass

class ConstantFunction(RBAFunction):
    name = 'constant'
    
    def __init__(self, parameters, variable = None):
        super(ConstantFunction, self).__init__(variable)
        self.CONSTANT = parameters['CONSTANT']

    def evaluate(self, x):
        return self.CONSTANT

class ExponentialFunction(RBAFunction):
    name = 'exponential'
    
    def __init__(self, parameters, variable = None):
        super(ExponentialFunction, self).__init__(variable)
        self.RATE = parameters['RATE']
        
    def evaluate(self, x):
        return numpy.exp(self.RATE*x)

class IndicatorFunction(RBAFunction):
    name = 'indicator'
    
    def __init__(self, parameters, variable = None):
        super(IndicatorFunction, self).__init__(variable)
        self.X_MIN = parameters['X_MIN']
        self.X_MAX = parameters['X_MAX']
        
    def evaluate(self, x):
        return (x > self.X_MIN) and (x < self.X_MAX)

class LinearFunction(RBAFunction):
    name = 'linear'
    
    def __init__(self, parameters, variable = None):
        super(LinearFunction, self).__init__(variable)
        self.X_MIN = parameters['X_MIN']
        self.X_MAX = parameters['X_MAX']
        self.LINEAR_COEF = parameters['LINEAR_COEF']
        self.LINEAR_CONSTANT = parameters['LINEAR_CONSTANT']
        self.Y_MIN = parameters['Y_MIN']
        self.Y_MAX = parameters['Y_MAX']
        
    def evaluate(self, x):
        x_eval = min(max(x, self.X_MIN), self.X_MAX)
        y = self.LINEAR_COEF * x_eval + self.LINEAR_CONSTANT
        return min(max(y, self.Y_MIN), self.Y_MAX)

class MichaelisMentenFunction(RBAFunction):
    name = 'michaelisMenten'
    
    def __init__(self, parameters, variable = None):
        super(MichaelisMentenFunction, self).__init__(variable)
        self.kmax = float(parameters['kmax'])
        self.Km = float(parameters['Km'])
        self.Y_MIN = None
        if 'Y_MIN' in parameters:
            self.Y_MIN = parameters['Y_MIN']
        
    def evaluate(self, x):
        y = self.kmax * x / (x + self.Km)
        if self.Y_MIN:
            return max(y, self.Y_MIN)
        else:
            return y

class MultiplicationFunction(RBAFunction):
    name = 'multiplication'
    
    def __init__(self, function_handles):
        super(MultiplicationFunction, self).__init__(None)
        self.OPERANDS = function_handles
        
    def evaluate(self, x):
        y = 1
        for op in self.OPERANDS:
            y *= op.evaluate(x)
        return y

# list of accepted function names and classes implementing them
simple = [ConstantFunction, LinearFunction, IndicatorFunction,
          ExponentialFunction, MichaelisMentenFunction]
aggregate = [MultiplicationFunction]
valid_fns = {c.name:c for c in simple}
valid_aggs = {c.name:c for c in aggregate}

def build_function(type_, params, variable = None):
    try:
        # retrieve class implementing function
        fn_class = valid_fns[type_]
    except KeyError:
        print('Unknown function type: ' + type_ + '. Valid types are: '
              + ', '.join(valid_fns.keys()))
        raise UserWarning('Invalid function.')
    try:
        return fn_class(params, variable)
    except KeyError as error:
        print('Missing parameter: ' + error.message)
        raise UserWarning('Invalid function.')

def build_aggregate(self, agg, known_functions):
    try:
        # retrieve class implementing aggregate
        agg_class = valid_aggs[agg.type]
    except KeyError:
        print('Unknown aggregate type: ' + agg.type + '. Valid types are: '
              + ', '.join(valid_aggs.keys()))
        raise UserWarning('Invalid aggregate.')
    try:
        # retrieve functions used in aggregate
        fn_handles = [known_functions[ref.function] \
                      for ref in agg.function_references]
    except KeyError as error:
        print('Unknown function: ' + error.args[0])
        raise UserWarning('Invalid aggregate.')
    return agg_class(fn_handles)
