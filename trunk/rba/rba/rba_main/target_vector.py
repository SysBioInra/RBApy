
import numpy

class TargetVector(object):
    def __init__(self, values, known_functions, def_value = None):
        self._values = []
        self._functions = []
        for v in values:
            try:
                self._values.append(float(v))
                self._functions.append(None)
            except (ValueError, TypeError):
                try:
                    self._functions.append(known_functions[v])
                    self._values.append(0)
                except KeyError:
                    if def_value is not None:
                        self._values.append(def_value)
                        self._functions.append(None)
                    else:
                        raise UserWarning('Invalid value: ' + v)
            
    def compute(self, mu):
        result = numpy.array(self._values, dtype='float')
        for i, fn in enumerate(self._functions):
            if fn is None: continue
            value = fn.evaluate(mu)
            result[i] = value
        return result
