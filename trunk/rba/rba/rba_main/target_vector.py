
import numpy

class TargetVector(object):
    def __init__(self, value_nodes, functions):
        self._values = []
        self._functions = []
        for n in value_nodes:
            self._values.append(n.value)
            self._functions.append([functions[f] for f in n.function_references])
            
    def compute(self, mu):
        result = numpy.array(self._values, dtype='float')
        for i, fn_list in enumerate(self._functions):
            if not(fn_list): continue
            value = 1
            for fn in fn_list: value *= fn.evaluate(mu)
            result[i] = value
        return result
