from __future__ import division

from numpy import array, cos, cross, dot, mat, pi, sin, sqrt, zeros
from numpy.linalg import norm
from ._math import rotz, rotzxz, asunits
from ._vector import state_vector
from .algorithm._geodetic import geodetic_latitude

class Frame(object):
    def __init__(object):
        pass


class InertialFrame(Frame):
    def __init__(self):
        pass
        
    def toinertial(self, stv, t):
        return stv
        
    def tolocal(self, stv, t):
        return stv


class GeocentricFrame(Frame):
    def __init__(self, inc, lonasc, argve):
        self._A = rotzxz(lonasc, inc, argve)
    
    def toinertial(self, stv, t):
        return state_vector(dot(self._A, stv.r).A1, dot(self._A, stv.v).A1)
    
    def tolocal(self, stv, t):
        return state_vector(dot(self._A.T, stv.r).A1, dot(self._A.T, stv.v).A1)


class GeocentricRotatingFrame(GeocentricFrame):
    def __init__(self, inc, lonasc, argve, sidereal_rate):
        GeocentricFrame.__init__(self, inc, lonasc, argve)
        self._w = sidereal_rate
    
    def toinertial(self, stv, t):
        A = dot(rotz(t*self._w), self._A)
        r = dot(A, stv.r).A1
        v = dot(A, stv.v + cross(array([0,0,self._w]), stv.r)).A1
        return state_vector(r, v)
    
    def tolocal(self, stv, t):
        A = dot(rotz(t*self._w), self._A).T
        r = dot(A, stv.r).A1
        v = dot(A, stv.v).A1 - cross(array([0,0,self._w]), r)
        return state_vector(r, v)


class PerifocalFrame(GeocentricFrame):
    def __init__(self, inc, lonasc, argpe):
        GeocentricFrame.__init__(self, inc, lonasc, argpe)
    
    def tostatevector(self, loc):
        rp = array(list(loc.r) + [0])
        vp = array(list(loc.v) + [0])
        return state_vector(dot(self._A, rp).A1, dot(self._A, vp).A1)
    
    def tolocalvector(self, stv):
        return dot(self._A.T, stv.r).A1[0:2], dot(self._A.T, stv.v).A1[0:2]


class OrbitalFrame(Frame):
    def __init__(self, orbit):
        self.orbit = orbit
    
    def toinertial(self, stv, t):
        return stv + self.orbit.statevector_by_time(t)
    
    def tolocal(self, stv, t):
        return stv - self.orbit.statevector_by_time(t)


class GeodeticFrame(GeocentricRotatingFrame):
    def __init__(self, Rp, Re, inc, lonasc, argve, sidereal_rate):
        GeocentricRotatingFrame.__init__(self, inc, lonasc, argve, sidereal_rate)
        self.f = (Rp + Re)/Re
        self.e = sqrt(2*self.f-self.f**2)
        self.Rp = Rp
        self.Re = Re
    
    def surface_height(self, lat):
        return self.Re/sqrt(1-self.e**2*sin(lat)**2)
    
    def uniti(self, lat, lon, t):
        st = self._w * t
        return array([-sin(st), cos(st), 0])
    
    def unitj(self, lat, lon, t):
        st = self._w * t
        return array([-sin(lat)*cos(st), -sin(lat)*sin(st), cos(lat)])
    
    def unitk(self, lat, lon, t):
        st = self._w * t
        return array([cos(lat)*cos(st), cos(lat)*sin(st), sin(lat)])
    
    def surface_inertial_vector(self, lat, lon, altitude, t, v=zeros(3)):
        Rlat = self.surface_height(lat)
        r0 = array([Rlat*cos(lat), 0, Rlat*(1-self.f)**2*sin(lat)]) + self.unitk(lat, 0, 0)
        v0 = cross(array([0,0,self._w]), r0) + v
        A = dot(rotz(t*self._w), self._A)
        return state_vector(dot(A.T,r0).A1, dot(A.T,v0).A1)
    
    def geodetic_llav(self, stv, t):
        lat, alt = geodetic_latitude(stv.r, self.Re, self.e)
        if stv[1] > 0:
            lon = arccos(stv.r[0]/norm(stv.r)) - ((t*self._w) % (2*pi))
        else:
            lon = 2*pi - arccos(stv.r[0]/norm(stv.r)) - ((t*self._w) % (2*pi))
        while lon < 0:
            lon += 2*pi
        A = mat([self.uniti(lat, lon, t), self.unitj(lat, lon, t), self.unitk(lat, lon, t)])
        v = dot(A,stv.v).A1
        return lat, lon, alt, v
    
    def altitude(self, stv):
        _, alt = geodetic_latitude(stv.r, self.Re, self.e)
        return alt


def inertial_frame():
    return InertialFrame()

def geodetic_frame(Rp, Re, inc, lonasc, argve, period):
    return GeodeticFrame(Rp,Re,inc,lonasc,argve,(2*pi)/period)

def parse_geodetic_frame(geodetic_expr):
    expr_list = [e.strip() for e in geodetic_expr[1:-1].split(',')]
    Re, Rp, inc, lonasc, argve, period = asunits(expr_list, ['m','m','rad','rad','rad','sec'])
    return GeodeticFrame(Rp, Re, inc, lonasc, argve, (2*pi)/period)

def geocentric_frame(inc, lonasc, argve):
    return GeocentricFrame(inc, lonasc, argve)

def perifocal_frame(inc, lonasc, argpe):
    return PerifocalFrame(inc, lonasc, argpe)

def orbital_frame(orbit):
    return OrbitalFrame(orbit)

def parse_orbital_frame(kepler_expr, u, epoch):
    from ._kepler import parse_kepler
    return OrbitalFrame(parse_kepler(kepler_expr, u, epoch))

