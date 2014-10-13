from __future__ import division

from numpy import array, arccos, arcsin, cos, cross, dot, mat, pi, sin, sqrt, zeros
from numpy.linalg import norm
from ._math import rotz, rotzxz, asunits, uniti, unitk, unitj
from ._persistant import PersistantObject
from ._vector import statevector
from .algorithm._geodetic import geodetic_latitude


class Frame(PersistantObject):
    def __init__(self):
        pass
    
    def uniti(self, t):
        raise NotImplementedError
    
    def unitj(self, t):
        raise NotImplementedError
    
    def unitk(self, t):
        raise NotImplementedError


class InertialFrame(Frame):
    def __init__(self):
        Frame.__init__(self)
        
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
        InertialFrame.__init__(self)
        self.origin = origin
    
    def toinertial(self, stv, t):
        return stv + self.origin
    
    def tolocal(self, stv, t):
        return stv - self.origin


class FunctionalDisplacementFrame(InertialFrame):
    def __init__(self):
        InertialFrame.__init__(self)
    
    def toinertial(self, stv, t):
        return stv + self._get_origin(t)
    
    def tolocal(self, stv, t):
        return stv - self._get_origin(t)
    
    def _get_origin(self, t):
        raise NotImplementedError


class ConstantOrientationFrame(Frame):
    def __init__(self, R):
        Frame.__init__(self)
        self._A = R.I
        self._AI = R
        self.mapvar('_A', 'A')
        self.mapvar('_AI', 'AI')
    
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
        Frame.__init__(self)
        self._A = A
        self._AI = A.I
        self._w = w
        self.mapvar('_A', 'A')
        self.mapvar('_AI', 'AI')
        self.mapvar('_w', 'w')
    
    def _tolocal_A(self, t):
        return dot(rotz(t*self._w), self._A)
    
    def toinertial(self, stv, t):
        A = self._tolocal_A(t)
        r = dot(A, stv.r).A1
        v = dot(A, stv.v + cross(unitk*self._w, stv.r)).A1
        return statevector(r, v)
    
    def tolocal(self, stv, t):
        A = self._tolocal_A(t).I
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
        ConstantOrientationFrame.__init__(self, rotzxz(argve, inc, lonasc))


class GeocentricRotatingFrame(ConstantRotationFrame):
    def __init__(self, inc, lonasc, argve, sidereal_rate):
        ConstantRotationFrame.__init__(self, rotzxz(argve, inc, lonasc), sidereal_rate)


class PerifocalFrame(GeocentricFrame):
    def __init__(self, inc, lonasc, argpe):
        GeocentricFrame.__init__(self, inc, lonasc, argpe)
    
    def tostatevector(self, loc):
        rp = array(list(loc.r) + [0])
        vp = array(list(loc.v) + [0])
        return statevector(rp, vp)
    
    def tolocalvector(self, stv):
        return stv.r[0:2], stv.v[0:2]


class OrbitalFrame(FunctionalDisplacementFrame):
    def __init__(self, orbit):
        FunctionalDisplacementFrame.__init__(self)
        self.orbit = orbit
    
    def _get_origin(self, t):
        return self.orbit.statevector_by_time(t)


class RotatingEllipsoide(PersistantObject):
    def __init__(self, Rp, Re, inc, lonasc, argve, sidereal_rate):
        self.frame = GeocentricRotatingFrame(inc, lonasc, argve, sidereal_rate)
        self.f = (Rp + Re)/Re
        self.e = sqrt(2*self.f-self.f**2)
        self.Rp = Rp
        self.Re = Re
        self._w = self.frame._w
        self.mapvar('_w','w')
    
    def surface_height(self, lat):
        return self.Re/sqrt(1-self.e**2*sin(lat)**2)
    
    def uniti(self, lat, lon, t):
        st = t * self._w + lon
        u = array([-sin(st), cos(st), 0])
        return dot(self.frame._A, u).A1
    
    def unitj(self, lat, lon, t):
        st = t * self._w + lon
        u = array([-sin(lat)*cos(st), -sin(lat)*sin(st), cos(lat)])
        return dot(self.frame._A, u).A1
    
    def unitk(self, lat, lon, t):
        st = t * self._w + lon
        u = array([cos(lat)*cos(st), cos(lat)*sin(st), sin(lat)])
        return dot(self.frame._A, u).A1
    
    def surface_inertial_statevector(self, lat, lon, altitude, t):
        rlat = self.surface_height(lat)
        r0 = rlat * array([cos(lat), 0, (1-self.f)**2*sin(lat)]) + altitude*self.unitk(lat, 0, 0)
        v0 = zeros(3)
        return self.frame.toinertial(rotz(lon) * statevector(r0, v0), t)
    
    def geodetic_llav(self, stv, t):
        lat, alt = geodetic_latitude(stv.r, self.Re, self.e)
        r, v = self.frame.tolocal(stv, t).rv
        d = norm(stv.r[0:2])
        x, y, _ = r/d
        if y >= 0:
            lon = arccos(x)
        else:
            lon = (2*pi - arccos(x)) % (2*pi)
        i, j, k = self.uniti(lat, lon, 0), self.unitj(lat, lon, 0), self.unitk(lat, lon, 0)
        A = mat([i, j, k])
        return lat, lon, alt, dot(A, v).A1
    
    #def altitude(self, stv):
    #    _, alt = geodetic_latitude(stv.r, self.Re, self.e)
    #    return alt


class SurfaceFrame(Frame):
    def __init__(self, surface, lat, lon, alt):
        self.surface = surface
        self.lat, self.lon, self.alt = lat, lon, alt
    
    def uniti(self, t):
        return self.surface.uniti(self.lat, self.lon, t)
    
    def unitj(self, t):
        return self.surface.unitj(self.lat, self.lon, t)
    
    def unitk(self, t):
        return self.surface.unitk(self.lat, self.lon, t)
    
    def tolocal(self, stv, t):
        origin = self.surface.surface_inertial_statevector(self.lat, self.lon, self.alt, t)
        A = mat([self.uniti(t), self.unitj(t), self.unitk(t)])
        return A*(stv - origin)
    
    def toinertial(self, stv, t):
        origin = self.surface.surface_inertial_statevector(self.lat, self.lon, self.alt, t)
        A = mat([self.uniti(t), self.unitj(t), self.unitk(t)])
        return A.I*stv + origin


def inertial_frame():
    return InertialFrame()


def surface_frame(surface, lat, lon, alt):
    return SurfaceFrame(surface, lat, lon, alt)


def parse_surface_frame(surface, lla_expr):
    expr_list = [e.strip() for e in lla_expr[1:-1].split(',')]
    lat, lon, alt = asunits(expr_list, ['rad', 'rad', 'm'])
    return surface_frame(surface, lat, lon, alt)


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

