from __future__ import division

from numpy import array, dot
from ._math import rotzxz
from ._vector import state_vector

class Frame(object):
    def __init__(object):
        pass


class PerifocalFrame(Frame):
    def __init__(self, inc, lonasc, argpe):
        self._A = rotzxz(lonasc, inc, argpe)
        self.inc = inc
        self.lonasc = lonasc
        self.argpe = argpe
    
    def tostatevector(self, loc):
        rp = array(list(loc.r) + [0])
        vp = array(list(loc.v) + [0])
        return state_vector(dot(self._A, rp).A1, dot(self._A, vp).A1)


class OrbitalFrame(Frame):
    def __init__(self, orbit):
        self.orbit = orbit
    
    def toinertial(self, stv, t):
        return stv + self.orbit.statevector_by_time(t)
    
    def tolocal(self, stv, t):
        return stv - self.orbit.statevector_by_time(t)


def perifocal_frame(inc, lonasc, argpe):
    return PerifocalFrame(inc, lonasc, argpe)

def orbital_frame(orbit):
    return OrbitalFrame(orbit)

