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
            self.M0 = M
        elif a < 0:
            self.mean_motion = sqrt(u/(-self.semi_major_axis)**3)
            self.M0 = M
     
    @classmethod
    def from_planet_paremters(cls, u, a, e, i, arg_pe, lon_asc, M, epoch=0.):
        return cls(u, a, e, i, lon_asc, arg_pe, M, epoch)
    
    @classmethod
    def from_rvu(cls, r, v, u, epoch=0.):
        evec = (1./u)*((norm(v)**2-u/norm(r))*r-dot(r,v)*v)
        e = norm(evec)
        a = 1/(2/norm(r)-norm(v)**2/u)
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
            if dot(r, v) < 0:
                M = 2*pi - M
        elif a < 0:
            F = arccosh((e+costa)/(1+e*costa))
            M = e*sinh(F) - F
            if dot(r, v) < 0:
                M = -M
        
        return cls(u, a, e, i, la, argpe, M, epoch)
    
    @classmethod
    def circular_equitorial(cls, u, r, lon_epoch, epoch=0.):
        return cls(u, r, 0.0, 0.0, lon_epoch, 0.0, 0.0, epoch)
        

    def mean_anomaly(self, t):
        return self.M0 + self.mean_motion * (t - self.epoch)
    
    def eccentric_anomaly(self, t):
        M = self.mean_anomaly(t)
        e = self.eccentricity
        f = lambda nE: nE - e*sin(nE) - M
        fp = lambda nE: 1 - e*cos(nE)
        fp2 = lambda nE: e*sin(nE)
        return newton(f,0,fp,fprime2=fp2)
    
    def hyperbolic_eccentric_anomaly(self, t):
        M = self.mean_anomaly(t)
        e = self.eccentricity
        f = lambda nF: e*sinh(nF) - nF - M
        fp = lambda nF: e*cosh(nF) - 1
        fp2 = lambda nF: -e*sinh(nF)
        return newton(f,0,fp,fprime2=fp2,maxiter=1000)
    
    def true_anomaly(self, t):
        e = self.eccentricity
        a = self.semi_major_axis
        if a > 0:
            E = self.eccentric_anomaly(t)
            ta = 2*arctan(sqrt((1.+e)/(1.-e))*tan(E/2.))
            return ta
        elif a < 0 and e != 1:
            F = self.hyperbolic_eccentric_anomaly(t)
            ta = arccos((e-cosh(F))/(e*cosh(F)-1))
            if F < 0:
                return -ta
            else:
                return ta
        raise Exception("Peribolic orbits are not supported")
    
    def prograde(self, t, theta=None):
        e = self.eccentricity
        a = self.semi_major_axis
        if theta is None:
            theta = self.true_anomaly(t)
        fa = 0.5*pi + theta - arctan(e*sin(theta)/(1+e*cos(theta)))
        return self.orientation * array([cos(fa), sin(fa), 0])
    
    def radialin(self, t, theta=None):
        e = self.eccentricity
        if theta is None:
            theta = self.true_anomaly(t)
        ct, st, _ = cossin(theta)
        rad = Ax(rotz(self.flightvec(t, theta)), -uniti)
        return self.orientation * rad
    
    def normal(self):
        return self.orientation * unitk
    
    def pe(self):
        p = self.semi_latus_rectum
        e = self.eccentricity
        return self.orientation * (uniti*p/(1+e))
    
    def ap(self):
        p = self.semi_latus_rectum
        e = self.eccentricity
        return self.orientation * (-uniti*p/(1-e))
    
    def time_to_pe(self, t):
        M = 2*pi - self.mean_anomaly(t)
        return M / self.mean_motion
    
    def time_to_ap(self, t):
        M = self.mean_anomaly(t)
        if M > pi:
            return self.time_to_pe(t) + (2*pi-M)/self.mean_motion
        else:
            return (pi-M) / self.mean_motion
    
    def r(self, t, theta=None):
        from .util import cossin
        from numpy import cos
        p,e = self.semi_latus_rectum, self.eccentricity
        if theta is None:
            theta = self.true_anomaly(t)
        l = p/(1+e*cos(theta))
        return self.orientation * ((cossin(theta))*l)
    
    def v(self, t, theta=None):
        p,e,a,u = self.semi_latus_rectum, self.eccentricity, self.semi_major_axis, self.GM
        if theta is None:
            theta = self.true_anomaly(t)
        ct, st, _ = cossin(theta)
        l = p/(1+e*ct)
        vel = sqrt(u*(2./l-1./a))
        return vel*self.prograde(t, theta)
    
    def rv(self, t, theta=None):
        if theta is None:
            theta = self.true_anomaly(t)
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
        
