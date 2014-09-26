from __future__ import division

from numpy import array, arcsin, cos, arccos, pi
from numpy.linalg import norm
from scipy.integrate import odeint, ode


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
    
    def _get_r(self):
        return self._vector[0]
    
    def _set_r(self, r):
        self._vector[0] = r
    
    def _get_v(self):
        return self._vector[1]
    
    def _set_v(self, v):
        self._vector[1] = v
    
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
    dar = property(_get_dar)


class PerifocalVector(RVVector):
    def __init__(self, r, v):
        RVVector.__init__(r,v)


def statevector(r,v):
    return StateVector(r,v)


def perifocal_vector(r,v):
    return PerifocalVector(r,v)


