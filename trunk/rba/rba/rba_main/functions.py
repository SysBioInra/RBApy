#
# File containing valid RBA functions.
#
import numpy
import abc

class Functions:
    def __init__(self, data):
        self.functions = {}
        for fn in data.parameters.functions:
            params = {p.id: p.value for p in fn.parameters}
            self.functions[fn.id] = build_function(fn.type, params, fn.variable)

    def __getitem__(self, fn_id):
        return self.functions[fn_id]

class RBAFunction(object):
    __metaclass__ = abc.ABCMeta

    def __init__(self, type_, parameters, variable = None):
        self.type_ = type_
        self.parameters = parameters
        self.variable = variable

    @abc.abstractmethod
    def evaluate(self, x): pass

class ConstantFunction(RBAFunction):
    def __init__(self, parameters, variable = None):
        super(ConstantFunction, self).__init__('constant', parameters, variable)
        self.CONSTANT = parameters['CONSTANT']

    def evaluate(self, x):
        return self.CONSTANT

class ExponentialFunction(RBAFunction):
    def __init__(self, parameters, variable = None):
        super(ExponentialFunction, self).__init__('exponential', parameters, variable)
        self.RATE = parameters['RATE']
        
    def evaluate(self, x):
        return numpy.exp(self.RATE*x)

class IndicatorFunction(RBAFunction):
    def __init__(self, parameters, variable = None):
        super(IndicatorFunction, self).__init__('indicator', parameters, variable)
        self.X_MIN = parameters['X_MIN']
        self.X_MAX = parameters['X_MAX']
        
    def evaluate(self, x):
        return (x > self.X_MIN) and (x < self.X_MAX)

class LinearFunction(RBAFunction):
    def __init__(self, parameters, variable = None):
        super(LinearFunction, self).__init__('linear', parameters, variable)
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
    def __init__(self, parameters, variable = None):
        super(MichaelisMentenFunction, self).__init__('michaelisMenten', parameters, variable)
        self.kmax = float(parameters['kmax'])
        self.Km = float(parameters['Km'])
        self.Y_MIN = None
        if parameters.has_key('Y_MIN'):
            self.Y_MIN = parameters['Y_MIN']
        
    def evaluate(self, x):
        y = self.kmax * x / (x + self.Km)
        if self.Y_MIN:
            return max(y, self.Y_MIN)
        else:
            return y

class MultiplicationFunction(RBAFunction):
    def __init__(self, function_ids, function_handles, variable = None):
        parameters = {'OPERANDS': function_ids}
        super(MultiplicationFunction, self).__init__('multiplication', parameters, variable)
        self.OPERANDS = function_handles
        
    def evaluate(self, x):
        y = 1
        for op in self.OPERANDS:
            y *= op.evaluate(x)
        return y

# list of accepted class types and classes implementing them
valid_types = {'constant': ConstantFunction,
               'linear': LinearFunction,
               'indicator': IndicatorFunction,
               'exponential': ExponentialFunction,
               'michaelisMenten': MichaelisMentenFunction}

def build_function(type_, parameters, variable = None):
    try:
        fn_class = valid_types[type_]
    except KeyError:
        print('Unknown function type: ' + type_ + '. Valid types are: '
              + ', '.join(valid_types.keys()))
        raise UserWarning('Invalid function.')
    try:
        return fn_class(parameters, variable)
    except KeyError as error:
        print('Missing parameter: ' + error.message)
        raise UserWarning('Invalid function.')
