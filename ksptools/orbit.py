from __future__ import division

class KeplerOrbit(object):
    def __init__(self, u, a, e, inc, lon_asc, arg_pe, M, epoch=0.):
        from numpy import sqrt, pi
        from .util import EulerAngle
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
    def from_rvu(cls, rvec, velvec, u, epoch=0.):
        from numpy import cross, array, dot, arccos, arccosh, sin, sinh, cos, pi
        from numpy.linalg import norm
        from .util import veccos
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

    def E_anom(self, t):
        from numpy import cos, sin
        from scipy.optimize import newton
        t -= self.epoch
        M = self.mean_anomally + (self.mean_motion * t)
        e = self.eccentricity
        f = lambda nE: nE - e*sin(nE) - M
        fp = lambda nE: 1 - e*cos(nE)
        fp2 = lambda nE: e*sin(nE)
        return newton(f,0,fp,fprime2=fp2)
    
    def F_anom(self, t):
        from numpy import cosh, sinh
        from scipy.optimize import newton
        t -= self.epoch
        M = self.mean_anomally + (self.mean_motion * dt)
        e = self.eccentricity
        f = lambda nF: e*sinh(nF) - nF - M
        fp = lambda nF: e*cosh(nF) - 1
        fp2 = lambda nF: e*sinh(nF)
        return newton(f,0,fp,fprime2=fp2)
    
    def true_anom(self, t):
        from numpy import arctan, arccos, tan, sqrt, cosh
        e = self.eccentricity
        a = self.semi_major_axis
        if a > 0:
            E = self.eccentricity_anom(t)
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
        from .util import cossin
        from numpy import array, arctan, arctan2, pi
        e = self.eccentricity
        if theta is None:
            theta = self.true_anom(t)
        ct, st, _ = cossin(theta)
        return array([-st, e+ct, 0])
        #return 0.
    
    def prograde(self, t, theta=None):
        from .util import cossin, Ax, rotz
        from numpy import array
        e = self.eccentricity
        if theta is None:
            theta = self.true_anom(t)
        ct, st, _ = cossin(theta)
        
        pgd = self.flightvec(t, theta)
        return self.orientation * pgd
    
    def radialin(self, t, theta=None):
        from .util import cossin, Ax, rotz
        from numpy import array
        e = self.eccentricity
        if theta is None:
            theta = self.true_anom(t)
        ct, st, _ = cossin(theta)
        rad = Ax(rotz(self.flightvec(t, theta)), array([-1., 0., 0.]))
        return self.orientation * rad
    
    def normal(self):
        from numpy import array
        return self.orientation * array([0.,0.,1.])
    
    def r(self, t, theta=None):
        from .util import cossin
        from numpy import cos
        p,e = self.semi_latus_rectum, self.eccentricity
        if theta is None:
            theta = self.true_anom(t)
        l = p/(1+e*cos(theta))
        return self.orientation * ((cossin(theta))*l)
    
    def pe(self):
        from numpy import array
        p = self.semi_latus_rectum
        e = self.eccentricity
        return self.orientation * ((array([1.,0.,0.]))*p/(1+e))
    
    def v(self, t, theta=None):
        from .util import cossin, Ax, rotz
        from numpy import array, arctan2, sqrt
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
        from numpy import pi
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
        
