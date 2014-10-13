from __future__ import division

from ._math import unit

from numpy import array, arcsin, cos, dot, arccos, pi, zeros
from numpy.linalg import norm
from scipy.interpolate import interp1d


class _RVVector(object):
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
    
    def _get_rv(self):
        return tuple(self._vector)
    
    @classmethod
    def _dims(cls):
        raise NotImplementedError
    
    @classmethod
    def zero(cls):
        return cls(zeros(cls._dims()), zeros(cls._dims()))
    
    @classmethod
    def interp(cls, ts, stvs):
        dims = list(range(cls._dims()))
        rifunc = [interp1d(ts, [s.r[i] for s in stvs], 'quadratic') for i in dims]
        vifunc = [interp1d(ts, [s.v[i] for s in stvs], 'quadratic') for i in dims]
        def interpfunc(t):
            return cls(
                array([rifunc[i](t) for i in dims]),
                array([vifunc[i](t) for i in dims]))
        return interpfunc
    
    def __str__(self):
        return '[{}*<{}>, {}*<{}>]'.format(
            norm(self.r), ','.join(map(str,unit(self.r))),
            norm(self.v), ','.join(map(str,unit(self.v))))
    
    def __getinitargs__(self):
        return tuple(self._vector)
    
    r = property(_get_r, _set_r)
    v = property(_get_v, _set_v)
    rv = property(_get_rv)


class _StateVector(_RVVector):
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


class _PerifocalVector(_RVVector):
    @classmethod
    def _dims(cls):
        return 2


statevector = _StateVector
perifocal_vector = _PerifocalVector


