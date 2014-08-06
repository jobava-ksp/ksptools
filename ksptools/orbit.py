from __future__ import division

from numpy import array, dot, cross, sqrt, pi
from numpy import cos, sin, tan, arcsin, arccos, arctan
from numpy import cosh, sinh, tanh, arccosh, arcsinh, arctanh
from numpy.linalg import norm
from scipy.optimize import newton
from .util import EulerAngle, veccos, vecsin, cossin, Ax, rotz, uniti, unitj, unitk


class KeplerOrbit(object):
    def __init__(self, u, a, e, inc, lon_asc, arg_pe, M, epoch=0.):
        self.GM = u
        self.eccentricity = e
        self.semi_latus_rectum = a*(1-e**2)
        self.semi_major_axis = a
        self.epoch = epoch
        self.orientation = EulerAngle.from_pts(lon_asc, inc, arg_pe)
        if a > 0:
            self.mean_motion = sqrt(u/self.semi_major_axis**3)
            self.mean_anomally = M
        elif a < 0:
            self.mean_motion = sqrt(u/(-self.semi_major_axis)**3)
            self.mean_anomally = M
     
    @classmethod
    def from_planet_paremters(cls, u, a, e, i, arg_pe, lon_asc, M, epoch=0.):
        return cls(u, a, e, i, lon_asc, arg_pe, M, epoch)
    
    @classmethod
    def from_rvu(cls, r, v, u, epoch=0.):
        evec = (1./u)*((dot(v,v)-u/norm(r))*r-dot(r,v)*v)
        e = norm(evec)
        a = 1/(2/norm(r)-dot(v,v)/u)
        p = a*(1-e**2)
        
        h = cross(r, v)
        i = arccos(h[2]/norm(h))
         
        n = cross(unitk, h)
        if norm(n) == 0:
            n = evec
        
        la = arccos(n[0]/norm(n))
        if n[1] < 0: la = 2*pi - la
        
        argpe = arccos(veccos(n,evec))
        if evec[2] < 0: argpe = 2*pi - argpe
        
        costa = veccos(evec, r)
        if a > 0:
            E = arccos((e+costa)/(1+e*costa))
            M = E - e*sin(E)
        elif a < 0:
            F = arccosh((e+costa)/(1+e*costa))
            M = e*sinh(F) - F
        if dot(r, v) < 0: M = 2*pi - M
        
        return cls(u, a, e, i, la, argpe, M, epoch)

    def E_anom(self, t):
        t -= self.epoch
        M = self.mean_anomally + (self.mean_motion * t)
        e = self.eccentricity
        f = lambda nE: nE - e*sin(nE) - M
        fp = lambda nE: 1 - e*cos(nE)
        fp2 = lambda nE: e*sin(nE)
        return newton(f,0,fp,fprime2=fp2)
    
    def F_anom(self, t):
        t -= self.epoch
        M = self.mean_anomally + (self.mean_motion * t)
        e = self.eccentricity
        f = lambda nF: e*sinh(nF) - nF - M
        fp = lambda nF: e*cosh(nF) - 1
        fp2 = lambda nF: e*sinh(nF)
        return newton(f,0,fp,fprime2=fp2)
    
    def true_anom(self, t):
        e = self.eccentricity
        a = self.semi_major_axis
        if a > 0:
            E = self.E_anom(t)
            return 2*arctan(sqrt((1.+e)/(1.-e))*tan(E/2.))
        elif a < 0 and e != 1:
            F = self.F_anom(t)
            ta = arccos((e-cosh(F))/(e*cosh(F)-1))
            if F < 0:
                return -ta
            else:
                return ta
        raise Exception("Peribolic orbits are not supported")
    
    def prograde(self, t, theta=None):
        e = self.eccentricity
        if theta is None:
            theta = self.true_anom(t)
        ct, st, _ = cossin(theta)
        return self.orientation * array([-st, e+ct, 0])
    
    def radialin(self, t, theta=None):
        e = self.eccentricity
        if theta is None:
            theta = self.true_anom(t)
        ct, st, _ = cossin(theta)
        rad = Ax(rotz(self.flightvec(t, theta)), -uniti)
        return self.orientation * rad
    
    def normal(self):
        return self.orientation * unitk
    
    def r(self, t, theta=None):
        from .util import cossin
        from numpy import cos
        p,e = self.semi_latus_rectum, self.eccentricity
        if theta is None:
            theta = self.true_anom(t)
        l = p/(1+e*cos(theta))
        return self.orientation * ((cossin(theta))*l)
    
    def pe(self):
        p = self.semi_latus_rectum
        e = self.eccentricity
        return self.orientation * (uniti*p/(1+e))
    
    def v(self, t, theta=None):
        p,e,a,u = self.semi_latus_rectum, self.eccentricity, self.semi_major_axis, self.GM
        if theta is None:
            theta = self.true_anom(t)
        ct, st, _ = cossin(theta)
        l = p/(1+e*ct)
        vel = sqrt(u*(2./l-1./a))
        return vel*self.prograde(t, theta)
    
    def rv(self, t, theta=None):
        if theta is None:
            theta = self.true_anom(t)
        return self.r(t, theta), self.v(t, theta)
    
    def period(self):
        if self.semi_major_axis > 0:
            return 2*pi/self.mean_motion
        elif self.semi_major_axis < 0:
            return float('inf')
        
    @classmethod
    def from_config(cls, conf_parser, section, u):
        from . import unit
        a = conf_parser.getfloat(section, 'a')
        e = conf_parser.getfloat(section, 'e')
        i = unit.tounit(conf_parser.get(section, 'i'), 'rad')
        arg_pe = unit.tounit(conf_parser.get(section, 'arg_pe'), 'rad')
        lon_asc = unit.tounit(conf_parser.get(section, 'lon_asc'), 'rad')
        M = unit.tounit(conf_parser.get(section, 'mean_anom'), 'rad')
        if conf_parser.has_option(section, 'epoch'):
            epoch = conf_parser.getfloat(section, 'epoch')
        else:
            epoch = 0.0
        return cls.from_planet_paremters(u, a, e, i, arg_pe, lon_asc, M, epoch)
        
