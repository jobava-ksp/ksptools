from __future__ import division

from numpy import array, arccos, arctan, cross, cos, dot, log, pi, sin, sinh, sqrt, tan, tanh

from numpy.linalg import norm
from scipy.optimize import newton

from ._math import C, S, asunits, uniti, unitk
from ._vector import perifocal_vector, statevector
from ._frame import perifocal_frame


def _herarp(stv, u):
    r, v = stv.r, stv.v
    h = norm(cross(r, v))
    r, v = nrom(r), norm(v)
    e = sqrt(1 + (h**2/u**2)*(v**2-2*u/r))
    rpe = (h**2/u)*(1/(1+e))
    rpa = (h**2/u)*(1/(1-e))
    return h, e, rpe, rap


class KeplerOrbit(object):
    def __init__(self, ecc, sma, inc, lonasc, argpe, u, epoch):
        """
        :type ecc: float
        :type sma: float
        :type inc: float
        :type lonasc: float
        :type argpe: float
        :type u: float
        :type epoch: float
        """
        p = sma*(1-ecc**2)
        h = sqrt(p*u)
        self.eccentricity = ecc
        self.semimajor_axis = sma
        self.semilatus_rectum = p
        self.specific_angular_momentum = h
        self.GM = u
        self.rpe = array([1,0,0]) * p/(1+ecc)
        self.rap = array([-1,0,0]) * p/(1-ecc)
        self.vpe = array([0,ecc + 1,0]) * (u/h)
        self.vap = array([0,ecc - 1,0]) * (u/h)
        self.frame = perifocal_frame(inc, lonasc, argpe)
        self.normal = self.frame.unitk(0)
        self.epoch = epoch
        if ecc < 1:
            self.period = 2*pi*sqrt(sma**3/u)
        else:
            self.period = float('+inf')
    
    @staticmethod
    def lambert(r0, t0, r1, t1, u):
        """
        :type r0: numpy.ndarray
        :type t0: float
        :type r1: numpy.ndarray
        :type t1: float
        :type u: float
        :rtype: KeplerOrbit
        """
        from .algorithm._lambert import lambert
        return lambert(r0, t0, r1, t1, u)
    
    @classmethod
    def from_statevector(cls, stv, u, epoch):
        """
        :type stv: ksptools._vector.StateVector:
        :type u: float
        :type epoch: float
        :rtype: KeplerOrbit
        """
        rad = norm(stv.r)
        vel = norm(stv.v)
        vrad = (dot(stv.v, stv.r)/dot(stv.r, stv.r))*stv.r
        h = cross(stv.r, stv.v)
        
        inc = arccos(h[2]/norm(h))
        n = cross(unitk,h)
        if norm(n) == 0:
            n = uniti
        if n[1] >= 0:
            lonasc = arccos(n[0]/norm(n))
        else:
            lonasc = 2*pi - arccos(n[0]/norm(n))
        
        e = (1/u)*((vel**2-u/rad)*stv.r-rad*vrad*stv.v)
        ecc = norm(e)
        
        if e[2] >= 0:
            argpe = arccos(dot(n,e)/(ecc*norm(n)))
        else:
            argpe = 2*pi - arccos(dot(n,e)/(ecc*norm(n)))

        if vrad[0] >= 0:
            ta = arccos(dot(e,stv.r)/(ecc*rad))
        else:
            ta = 2*pi - arccos(dot(e,stv.r)/(ecc*rad))
        
        slr = dot(h,h)/u
        sma = slr/(1-ecc**2)
        if ecc == 1:
            mm = u**2/norm(h)**3
            Mp = 0.5*tan(ta/2) + (1/6)*tan(ta/2)**3
            dt = Mp/mm
        elif ecc < 1:
            E = 2*arctan(sqrt((1-ecc)/(1+ecc))*tan(ta/2))
            dt = (E - ecc*sin(E))/sqrt(u/sma**3)
        elif ecc > 1:
            F = log((sqrt(ecc+1)+sqrt(ecc-1)*tan(ta/2))/(sqrt(ecc+1)-sqrt(ecc-1)*tan(ta/2)))
            dt = (ecc*sinh(F)-F)/sqrt(u/(-sma)**3)
        return cls(ecc, sma, inc, lonasc, argpe, u, epoch-dt)
    
    #@classmethod
    #def from_rvu(cls, r, v, u, epoch):
    #   return cls.from_statevector(statevector(r,v), u, epoch)
    
    @classmethod
    def from_parameters(cls, sma, ecc, inc, argpe, lonasc, M0, u, epoch):
        if ecc == 1:
            raise NotImplementedError
        elif ecc < 1:
            dt = M0/sqrt(u/sma**3)
        elif ecc > 1:
            dt = M0/sqrt(u/(-sma)**3)
        return cls(ecc, sma, inc, lonasc, argpe, u, epoch-dt)
    
    def perifocal_by_ta(self, ta):
        e, p, h, u = self.eccentricity, self.semilatus_rectum, self.specific_angular_momentum, self.GM
        ru = array([cos(ta), sin(ta)])
        vu = array([-sin(ta), e + cos(ta)])
        return perifocal_vector((p/(1+e*cos(ta)))*ru, (u/h)*vu)
    
    def perifocal_by_time(self, time):
        return self.perifocal_by_ta(self.true_anomaly_by_time(time))
    
    def statevector_by_ta(self, ta):
        return self.frame.tostatevector(self.perifocal_by_ta(ta))
    
    def statevector_by_time(self, time):
        return self.frame.tostatevector(self.perifocal_by_time(time))
    
    def true_anomaly_by_time(self, time):
        X = self._universal_anomaly(time)
        e = self.eccentricity
        if e == 1:
            raise NotImplementedError
        elif e < 1:
            E = X/sqrt(self.semimajor_axis)
            ta = arccos((e - cos(E))/(e*cos(E)-1))
        elif e > 1:
            F = X/sqrt(-self.semimajor_axis)
            ta = 2*arctan(sqrt((e+1)/(e-1))*tanh(F/2))
        return ta
    
    def true_anomaly_by_vector(self, iv):
        r, v = self.frame.tolocalvector(iv)
        ct = dot(r,array([1,0]))
        if r[1] >= 0:
            return arccos(ct)
        else:
            return 2*pi - arccos(ct)
    
    def true_anomaly_by_distance(self, r):
        p, e = self.semilatus_rectum, self.eccentricity
        return arccos((p/r-1)/e)
    
    def time_by_ta(self, ta, ts):
        offset = self.time_anomaly_by_ta(ta) + self.epoch
        n = int((ts - offset)/self.period)
        return offset + (n+1)*self.period
    
    def _universal_anomaly(self, time):
        delt = time - self.epoch
        r0, vr0, u = self.rpe[0], self.vpe[1], self.GM
        a = 1/self.semimajor_axis
        
        def func(x):
            z = a*x**2
            return (r0*vr0/sqrt(u))*x**2*C(z) + (1-a*r0)*x**3*S(z) + r0*x - sqrt(u)*delt
        
        def gfunc(x):
            z = a*x**2
            return (r0*vr0/sqrt(u))*x*(1-a*x**2*S(z)) + (1-a*r0)*x**2*C(z) + r0
        
        return newton(func, abs(a)*sqrt(u)*delt, gfunc)
    
    def time_anomaly_by_ta(self, ta):
        e, a, u = self.eccentricity, self.semimajor_axis, self.GM
        if e == 1:
            raise NotImplementedError
        elif e < 1:
            E = 2*arctan(sqrt((1-e)/(1+e))*tan(ta/2))
            return (E - e*sin(E))/sqrt(u/a**3)
        elif e > 1:
            F = log((sqrt(e+1)+sqrt(e-1)*tan(ta/2))/(sqrt(e+1)-sqrt(e-1)*tan(ta/2)))
            return (e*sinh(F)-F)/sqrt(u/(-a)**3)
    
    #def _sqrtu_dt(self, x):
    #    r0, vr0, u = self.rpe[0], self.vpe[1], self.GM
    #    a = 1/self.semimajor_axis
    #    z = a*x**2
    #    return (r0*vr0/sqrt(u))*x**2*C(z) + (1-a*r0)*x**3*S(z) + r0*x
    scalar_herpra = staticmethod(_herpra)


def parse_kepler(kepler_expr, u, epoch):
    params = list(asunits(kepler_expr[1:-1].split(','),['m',None,'rad','rad','rad','rad'])) + [u, epoch]
    return KeplerOrbit.from_parameters(*params)


