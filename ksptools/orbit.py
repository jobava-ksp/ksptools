from __future__ import division

from numpy import array, dot, cross, sqrt, pi
from numpy import cos, sin, tan, arcsin, arccos, arctan
from numpy import cosh, sinh, tanh, arccosh, arcsinh, arctanh
from numpy.linalg import norm
from scipy.optimize import newton
from .util import EulerAngle, veccos, vecsin, cossin, Ax, rotz

zeropos = array([1.,0.,0.,])

def isnotNone(*args):
    return all(x is not None for x in args)

def spec_energy(u,r=None,a=None,e=None,v=None):
    # 1. u, v, r
    if isnotNone(u, v, r):
        return dot(v,v)/2 - u/norm(r)
    
    # 2. u, h, e
    if isnotNone(u, h, e):
        return -(u**2/(2*h**2))*(1-e**2)
    
    # 3. a, u
    if isnotNone(u, a):
        return -u/2*a
    
    return NotImplementedError

def spec_angular_momentum(r, v):
    return norm(cross(r, v))

def semi_major_axis(u,r=None,v=None,se=None):
    if isnotNone(u,r,v):
        return norm(r)*u/(2*u-norm(r)*dot(v,v))
    if isnotNone(u,sec):
        return -u/(2*sec)
        
    return NotImplementedError

def eccentricity(u, se=None, sam=None, r=None, v=None, a=None):
    # a. a, u
    if isnotNone(a,u):
        se = -u/2*a
    
    # 1. se, sam
    if isnotNone(se, sam):
        return sqrt(1+(2*se*sam**2)/u**2)
    
    # 2. se, r, v (r and v must be vectors)
    if isnotNone(se, r, v):
        return sqrt(1+(2*norm(cross(r, v))**2)/u**2)
    
    return NotImplementedError

def true_anomaly(a=None,e=None,r=None,E=None,F=None,dt=None,mm=None):
    # 1. a, e, r
    if isnotNone(a,e,r):
        if e != 1. and e > 0.0:
            ta = arccos((a*(1-e**2)-norm(r))/(norm(r)*e))
        elif e == 0.
            ta = veccos(r,array([1.,0.,0.]))
        if e < 1. and ta < 0.:
            ta = 2*pi - ta
        return ta
    
    # 2. e, E
    if isnotNone(e,E):
        return 2*arctan(sqrt((1.+e)/(1.-e))*tan(E/2.))
    
    # 3. e, F
    if isnotNone(e,F):
        ta = arccos((e-cosh(F))/(e*cosh(F)-1))
        if F < 0:
            return -ta
        else:
            return ta
    
    return NotImplementedError

def eccentric_anomaly(e, ta=None, dt=None, M0=None, mm=None, a=None, u=None, r=None):
    # a. a, u
    if isnotNone(a,u):
        if a > 0:
            mm = sqrt(a**3/u)
        elif a < 0:
            mm = sqrt((-a)**3/u)
        else:
            raise NotImplementedError
    
    # b. a, e, r
    if isnotNone(a,e,r):
        ta = true_anomaly(a, e, r)
    
    # 2. e, ta
    if isnotNone(e,ta):
        if e == 0:
            return ta
        elif e < 1:
            E = arccos((e+cos(ta))/(1+e*cos(ta)))
            if E < 0:
                E = 2*pi - E
            return E
        elif e > 1:
            F = arccosh((e+cos(ta))/(1+e*cos(ta)))
            return F
    
    # 3. e, dt, M0, mm
    if isnotNone(e,dt,M0,mm):
        M = M0 + mm*dt
        if e < 1:
            f = lambda nE: nE - e*sin(nE) - M
            fp = lambda nE: 1 - e*cos(nE)
            fp2 = lambda nE: e*sin(nE)
        elif e > 1:
            f = lambda nF: e*sinh(nF) - nF - M
            fp = lambda nF: e*cosh(nF) - 1
            fp2 = lambda nF: e*sinh(nF)
        else:
            raise NotImplementedError
        return newton(f,0,fp,fprime2=fp2)

    return NotImplementedError

def time_from_pe(u=None, a=None, e=None, EF=None, mm=None, ta=None, r=None):
    # a. a, e, r
    if isnotNone(a,e,r):
        ta = true_anomaly(a=a,r=r,e=e)
    
    # b. u, a
    if isnotNone(u,a):
        if a > 0:
            mm = sqrt(a**3/u)
        elif a < 0:
            mm = sqrt((-a)**3/u)
    
    # b. ta, e
    if isnotNone(ta, e):
        EF = eccentric_anomally(e=e, ta=ta)
    
    # 1. EF, e, mm
    if isnotNone(EF,e,m):
        if e < 1:
            return (EF - e*sin(EF))*mm
        elif e > 1:
            return (e*sinh(EF) - EF)*mm
    
    return NotImplementedError
    

class KeplerOrbit(object):
    def __init__(self, u, a, e, inc, lon_asc, arg_pe, M, epoch):
        self.u = u
        self.e = e
        self.p = a*(1-e**2)
        self.a = a
        self.orient = EulerAngle.from_pts(lon_asc, inc, arg_pe)
        self.epoch = epoch
        if a > 0:
            self.mean_motion = sqrt(u/self.a**3)
            self.mean_anom = M
        elif a < 0:
            self.mean_motion = sqrt(u/(-self.a)**3)
            self.mean_anom = M
     
    @classmethod
    def from_planet_paremters(cls, u, a, e, i, arg_pe, lon_asc, M, epoch):
        return cls(u, a, e, i, lon_asc, arg_pe, M, epoch)
    
    @classmethod
    def from_rvu(cls, rvec, velvec, u, epoch=0.):
        r = norm(rvec)
        v = norm(velvec)
        evec = (1./u)*((v**2-u/r)*rvec-dot(rvec,velvec)*velvec)
        e = norm(evec)
        a = 1/(2/r-v**2/u)
        p = a*(1-e**2)
        
        h = cross(rvec,velvec)
        i = arccos(h[2]/norm(h))
         
        n = cross(array([0.,0.,1.]), h)
        if norm(n) == 0:
            n = evec
        
        la = arccos(n[0]/norm(n))
        if n[1] < 0: la = 2*pi - la
        
        #_a = dot(n,evec)/(norm(n)*norm(evec))
        #_a = max(-1.0, min(1.0, _a))
        argpe = arccos(veccos(n,evec))
        if evec[2] < 0: argpe = 2*pi - argpe
        
        costa = veccos(evec, rvec)
        #costa = dot(evec,rvec)/(norm(evec)*norm(rvec))
        #_a = (e+costa)/(1+e*costa)
        #_a = max(-1.0, min(1.0, _a))
        if a > 0:
            E = arccos((e+costa)/(1+e*costa))
            M = E - e*sin(E)
        elif a < 0:
            F = arccosh((e+costa)/(1+e*costa))
            M = e*sinh(F) - F
        if dot(rvec, velvec) < 0: M = 2*pi - M
        
        return cls(u, a, e, i, la, argpe, M, epoch)

    def E_anom(self, dt):
        M = self.mean_anom + (self.mean_motion * dt)
        e = self.e
        f = lambda nE: nE - e*sin(nE) - M
        fp = lambda nE: 1 - e*cos(nE)
        fp2 = lambda nE: e*sin(nE)
        return newton(f,0,fp,fprime2=fp2)
    
    def F_anom(self, dt):
        M = self.mean_anom + (self.mean_motion * dt)
        e = self.e
        f = lambda nF: e*sinh(nF) - nF - M
        fp = lambda nF: e*cosh(nF) - 1
        fp2 = lambda nF: e*sinh(nF)
        return newton(f,0,fp,fprime2=fp2)
    
    def true_anom(self, t):
        t -= self.epoch
        e = self.e
        a = self.a
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
    
    def flightvec(self, t, theta=None):
        e = self.e
        if theta is None:
            theta = self.true_anom(t)
        ct, st, _ = cossin(theta)
        return array([-st, e+ct, 0])
    
    def prograde(self, t, theta=None):
        e = self.e
        if theta is None:
            theta = self.true_anom(t)
        ct, st, _ = cossin(theta)
        
        pgd = self.flightvec(t, theta)
        return self.orient.rotate(pgd)
    
    def radialin(self, t, theta=None):
        e = self.e
        if theta is None:
            theta = self.true_anom(t)
        ct, st, _ = cossin(theta)
        rad = Ax(rotz(self.flightvec(t, theta)), array([-1., 0., 0.]))
        return self.orient.rotate(rad)
    
    def normal(self):
        return self.orient.rotate(array([0.,0.,1.]))
    
    def r(self, t, theta=None):
        p,e = self.p, self.e
        if theta is None:
            theta = self.true_anom(t)
        l = p/(1+e*cos(theta))
        return self.orient.rotate(cossin(theta))*l
    
    def pe(self):
        p = self.p
        e = self.e
        return self.orient.rotate(array([1.,0.,0.]))*p/(1+e)
    
    def v(self, t, theta=None):
        p,e,a,u = self.p, self.e, self.a, self.u
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
        if self.a > 0:
            return 2*pi/self.mean_motion
        elif self.a < 0:
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
        
