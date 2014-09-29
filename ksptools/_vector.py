from __future__ import division

from ._math import unit

from numpy import array, arcsin, cos, dot, arccos, pi, zeros
from numpy.linalg import norm


class RVVector(object):
    def __init__(self, r, v):
        self._vector = array([r,v])
    
    def __add__(self, other):
        return type(self)(self.r + other.r, self.v + other.v)
    
    def __sub__(self, other):
        return type(self)(self.r - other.r, self.v - other.v)
    
    def __iadd__(self, other):
        self._vector += other._vector
    
    def __isub__(self, other):
        self._vector -= other._vector

    def __rmul__(self, other):
        return type(self)(dot(other, self.r).A1, dot(other, self.v).A1)

    def _get_r(self):
        return self._vector[0]
    
    def _set_r(self, r):
        self._vector[0] = r
    
    def _get_v(self):
        return self._vector[1]
    
    def _set_v(self, v):
        self._vector[1] = v
    
    @classmethod
    def _dims(cls):
        raise NotImplementedError
    
    @classmethod
    def zero(cls):
        return cls(zeros(cls._dims()), zeros(cls._dims()))
    
    def __str__(self):
        return '[{}*<{}>, {}*<{}>]'.format(
            norm(self.r), ','.join(map(str,unit(self.r))),
            norm(self.v), ','.join(map(str,unit(self.v))))
    
    r = property(_get_r, _set_r)
    v = property(_get_v, _set_v)


class StateVector(RVVector):
    def _get_dar(self):
        s = norm(self.r)
        l, m, n = self.r/s
        decl = arcsin(n)
        if m > 0:
            return decl, arccos(l/cos(decl)), s
        else:
            return decl, 2*pi - arccos(l/cos(decl)), s
    
    @classmethod
    def _dims(cls):
        return 3
    
    dar = property(_get_dar)


class PerifocalVector(RVVector):
    @classmethod
    def _dims(cls):
        return 2


statevector = StateVector
perifocal_vector = PerifocalVector


