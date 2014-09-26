from __future__ import division

from numpy import array, cos, cross, dot, mat, pi, sin, sqrt, zeros
from numpy.linalg import norm
from ._math import rotz, rotzxz, asunits, uniti, unitk, unitj
from ._vector import statevector
from .algorithm._geodetic import geodetic_latitude

class Frame(object):
    def __init__(object):
        pass
    
    def uniti(self, t):
        raise NotImplementedError
    
    def unitj(self, t):
        raise NotImplementedError
    
    def unitk(self, t):
        raise NotImplementedError


class InertialFrame(Frame):
    def __init__(self):
        pass
        
    def toinertial(self, stv, t):
        return stv
        
    def tolocal(self, stv, t):
        return stv
    
    def uniti(self, t):
        return uniti
    
    def unitj(self, t):
        return unitj
    
    def unitk(self, t):
        return unitk


class ConstantDisplacementFrame(InertialFrame):
    def __init__(self, origin):
        self.origin = origin
    
    def toinertial(self, stv, t):
        return stv + self.origin
    
    def tolocal(self, stv, t):
        return stv - self.origin


class FunctionalDisplacementFrame(InertialFrame):
    def __init__(self, origin_func):
        self._rvfunc = origin_func
    
    def toinertial(self, stv, t):
        origin = self._rvfunc(t)
        return stv + origin
    
    def tolocal(self, stv, t):
        origin = self._rvfunc(t)
        return stv - origin


class ConstantOrientationFrame(Frame):
    def __init__(self, A):
        self._A = A
        self._AI = A.I
    
    def toinertial(self, stv, t):
        return statevector(dot(self._AI, stv.r).A1, dot(self._AI, stv.v).A1)
    
    def tolocal(self, stv, t):
        return statevector(dot(self._A, stv.r).A1, dot(self._A, stv.v).A1)
    
    def uniti(self, t):
        return self._A[0].A1
    
    def unitj(self, t):
        return self._A[1].A1
    
    def unitk(self, t):
        return self._A[2].A1


class ConstantRotationFrame(Frame):
    def __init__(self, A, w):
        self._A = A
        self._AI = A.I
        self._w = w
    
    def _tolocal_A(self, t):
        return dot(rotz(-t*self._w), self._A)
    
    def toinertial(self, stv, t):
        A = self._tolocal_A(t).I
        r = dot(A, stv.r).A1
        v = dot(A, stv.v + cross(unitk*self._w, stv.r)).A1
        return statevector(r, v)
    
    def tolocal(self, stv, t):
        A = self._tolocal_A(t)
        r = dot(A, stv.r).A1
        v = dot(A, stv.v).A1 - cross(unitk*self._w, r)
        return statevector(r, v)
    
    def uniti(self, t):
        return self._tolocal_A(t)[0].A1
    
    def unitj(self, t):
        return self._tolocal_A(t)[1].A1
    
    def unitk(self, t):
        return self._tolocal_A(t)[2].A1


class GeocentricFrame(ConstantOrientationFrame):
    def __init__(self, inc, lonasc, argve):
        ConstantOrientationFrame.__init__(self, rotzxz(lonasc, inc, argve))


class GeocentricRotatingFrame(ConstantRotationFrame):
    def __init__(self, inc, lonasc, argve, sidereal_rate):
        ConstantRotationFrame.__init__(self, rotzxz(inc, lonasc, argve), sidereal_rate)


class PerifocalFrame(GeocentricFrame):
    def __init__(self, inc, lonasc, argpe):
        GeocentricFrame.__init__(self, inc, lonasc, argpe)
    
    def tostatevector(self, loc):
        rp = array(list(loc.r) + [0])
        vp = array(list(loc.v) + [0])
        return statevector(dot(self._A, rp).A1, dot(self._A, vp).A1)
    
    def tolocalvector(self, stv):
        return dot(self._A.T, stv.r).A1[0:2], dot(self._A.T, stv.v).A1[0:2]


class OrbitalFrame(FunctionalDisplacementFrame):
    def __init__(self, orbit):
        FunctionalDisplacementFrame.__init__(self, self._displacement)
        self.orbit = orbit
    
    def _displacement(self, t):
        return orbit._statevector_by_time(t)


class RotatingEllipsoide(object):
    def __init__(self, Rp, Re, inc, lonasc, argve, sidereal_rate):
        self.frame = GeocentricRotatingFrame(inc, lonasc, argve, sidereal_rate)
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
        return statevector(dot(A.T,r0).A1, dot(A.T,v0).A1)
    
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
    return RotatingEllipsoide(Rp,Re,inc,lonasc,argve,(2*pi)/period)

def parse_geodetic_frame(geodetic_expr):
    expr_list = [e.strip() for e in geodetic_expr[1:-1].split(',')]
    Re, Rp, inc, lonasc, argve, period = asunits(expr_list, ['m','m','rad','rad','rad','sec'])
    return RotatingEllipsoide(Rp, Re, inc, lonasc, argve, (2*pi)/period)

def geocentric_frame(inc, lonasc, argve):
    return GeocentricFrame(inc, lonasc, argve)

def perifocal_frame(inc, lonasc, argpe):
    return PerifocalFrame(inc, lonasc, argpe)

def orbital_frame(orbit):
    return OrbitalFrame(orbit)

def parse_orbital_frame(kepler_expr, u, epoch):
    from ._kepler import parse_kepler
    return OrbitalFrame(parse_kepler(kepler_expr, u, epoch))

