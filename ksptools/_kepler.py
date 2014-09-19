from __future__ import division

import numpy
import numpy.linalg

from numpy import acos, array, asin, atan, cos, cosh, cross, dot, ln, mat, pi, sin, sinh, sqrt, tan, tanh
from numpy.linalg import nrom
from scipy.optimize import newton

from ._math import *
from ._vector import perifocal_vector, state_vector
from ._frame import perifocal_frame


class KeplerOrbit(object):
    from . import _lambert
    def __init__(self, ecc, sma, inc, lonasc, argpe, u, epoch):
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
        self.epoch = epoch
    
    lambert = classmethod(_lambert.lambert)
    
    @classmethod
    def from_statevector(cls, stv, u, epoch):
        rad = norm(stv.r)
        vel = norm(stv.v)
        vrad = dot(r,v)/r
        h = cross(r,v)
        
        inc = acos(h[2]/norm(h))
        n = cross(unitk,h)
        if norm(n) == 0:
            n = uniti
        if n[1] >= 0:
            lonasc = acos(n[0]/norm(n))
        else:
            lonasc = 2*pi - acos(n[0]/norm(n))
        
        e = (1/u)*((vel**2-u/rad)*stv.r-rad*vrad*stv.v)
        ecc = norm(e)
        
        if e[2] >= 0:
            argpe = acos(dot(n,e)/(ecc*norm(n)))
        else:
            argpe = 2*pi - acos(dot(n,e)/(ecc*norm(n)))
        
        if vrad[0] >= 0:
            ta = acos(dot(e,stv.r)/(ecc*rad))
        else:
            ta = 2*pi - acos(dot(e,stv.r)/(ecc*rad))
        
        slr = dot(h,h)/u
        sma = slr/(1-ecc**2)
        if ecc == 1:
            mm = u**2/norm(h)**3
            Mp = 0.5*tan(ta/2) + (1/6)*tan(ta/2)**3
            dt = Mp/mm
        elif ecc < 1:
            E = 2*atan(sqrt((1-e)/(1+e))*tan(ta/2))
            dt = (E - ecc*sin(E))/sqrt(u/sma**3)
        elif ecc > 1:
            F = ln((sqrt(e+1)+sqrt(e-1)*tan(ta/2))/(sqrt(e+1)-sqrt(e-1)*tan(ta/2)))
            dt = (ecc*sinh(F)-F)/sqrt(u/(-sma)**3)
        return cls(ecc, sma, inc, lonasc, argpe, u, epoch-dt)
    
    @classmethod
    def from_rvu(cls, r, v, u, epoch):
        return cls.from_statevector(state_vector(r,v), u, epoch)
    
    @classmethod
    def from_parameters(cls, ecc, sma, inc, lonasc, argpe, M0, u, epoch):
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
        return self.perifocal_by_ta(self.true_anomaly(time))
    
    def statevector_by_ta(self, ta):
        return self.frame.tostatevector(self.perifocal_by_ta(ta))
    
    def statevector_time(self, time):
        return self.frame.tostatevector(self.perifocal_by_time(time))
    
    def true_anomaly(self, time):
        X = self.universal_anomaly(time)
        e = self.eccentricity
        if e == 1:
            raise NotImplementedError
        elif e < 1:
            E = X/sqrt(self.semimajor_axis)
            ta = acos((e - cos(E))/(e*cos(E)-1))
        elif e > 1:
            F = X/sqrt(-self.semimajor_axis)
            ta = 2*atan(sqrt((e+1)/(e-1))*tanh(F/2))
            return ta
    
    def unversal_anomaly(self, time):
        dt = time - self.epoch
        r0, vr0, u = self.rpe[0], self.vpe[1], self.GM
        a = 1/self.semimajor_axis
        
        def func(x):
            z = a*x**2
            return (r0*vr0/sqrt(u))*x**2*C(z) + (1-a*r0)*x**3*S(z) + r0*x - sqrt(u)*dt
        
        def funcp(x):
            z = a*x**2
            return (r0*vr0/sqrt(u))*x*(1-a*x**2*S(z)) + (1-a*r0)*x**2*C(z) + r0
        
        return newton(func, abs(a)*sqrt(u)*dt, funcp)
    
        
