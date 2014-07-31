from numpy import arccos, arcsin, arctan, arctan2, array, cross, cos, dot, mat, sin, sqrt, tan
from numpy.linalg import norm
from scipy.optimize import newton

from . import unit
from . import util

#TODO: move euler angle

class EulerAngle(object):    
    def __init__(self, p, t, s):
        Ap = util.rotz(p)  
        At = util.rotx(t)
        As = util.rotz(s)
        
        self.phi, self.theta, self.psi = p, t, s
        self.Ap, self.At, self.As = Ap, At, As
        self.A = Ap*At*As
    
    @classmethod
    def from_pts(cls,p,t,s):
        return cls(p,t,s)
    
    def rotate(self, r):
        return util.Ax(self.A,r)


class KeplerOrbit(object):
    def __init__(self, u, a, e, inc, lon_asc, arg_pe, M, epoch):
        self.u = u
        self.e = e
        self.a = a
        self.orient = EulerAngle.from_pts(lon_asc, inc, arg_pe)
        self.mean_motion = sqrt(u/self.a**3)
        self.mean_anom = (M - sqrt(u/a**3)*epoch) % pi
     
    @classmethod
    def from_planet_paremters(cls, u, a, e, i, arg_pe, lon_asc, M, epoch):
        return cls(u, a, e, i, lon_asc, arg_pe, M, epoch)
    
    #@classmethod
    #def from_rvu(cls, r, v, u, epoch):
    #    p,e,M,inc,lon_asc,arg_pe = KeplerOrbit._from_rvu(r,v,u)
    #    return cls(u, p, e, inc, lon_asc, arg_pe, M, epoch)

    @staticmethod
    def from_rvu(r, v, u, epoch=0.):
        vr = util.project(v,r)
        vt = util.reject(v,r)
        p = (norm(r)*norm(vt))**2/u
        e = sqrt(1 + 2*(norm(v)**2/2 - u/norm(r))*(norm(cross(r,v)))**2/u**2)
        if e == 0.0:
            theta = 0
        else:
            v0 = sqrt(u/p)
            acos_arg = min(1,max(-1,norm(vt)/(v0*e)-1/e))
            theta = arccos(acos_arg)
        E = arccos((e+cos(theta))/(1+e*cos(theta)))
        M = E - e*sin(E)
        
        x = cos(theta)*(r/norm(r)) - sin(theta)*(vt/norm(vt))
        y = sin(theta)*(r/norm(r)) + cos(theta)*(vt/norm(vt))
        z = cross(x,y)
        lon_asc = arctan2(z[0], -z[1])
        inc = arctan2(sqrt(z[0]**2+z[1]**2), z[2])
        arg_pe = arctan2(x[2],y[2])

        #return p,e,M,inc,lon_asc,arg_pe
        return cls(u, p, e, inc, lon_asc, arg_pe, M, epoch)

    def E_anom(self, dt):
        M = self.mean_anom + (self.mean_motion * dt)
        e = self.e
        n = self.mean_motion
        f = lambda nE: nE - e*sin(nE) - M
        fp = lambda nE: 1 - e*cos(nE)
        fp2 = lambda nE: e*sin(nE)
        return newton(f,0,fp,fprime2=fp2)
    
    def true_anom(self, t):
        E = self.E_anom(t)
        e = self.e
        return 2*arctan(sqrt((1+e)/(1-e))*tan(E/2))
    
    def flight_angle(self, t):
        from .util import cossin
        e = self.e
        theta = self.true_anom(t)
        ect, est, _ = e*cossin(t)
        return arctan2(1+ect, est)
    
    def prograde(self, t):
        from .util import cossin, Ax, rotz
        e = self.e
        theta = self.true_anom(t)
        ct, st, _ = cossin(theta)
        gamma = arctan2(1+e*ct,e*st)
        pgd = Ax(rotz(gamma + theta, array([-st, ct, 0])))
        return self.orient.rotate(pgd)
    
    def radialin(self, t):
        from .util import cossin
        theta = self.true_anom(t)
        return -self.orient.rotate(cossin(theta))
    
    def normal(self, t):
        return self.orient.rotate(array([0.,0.,1.]))
    
    def r(self, t):
        from .util import cossin
        p = self.p
        e = self.e
        theta = self.true_anom(t)
        l = p/(1+e*cos(theta))
        return self.orient.rotate(cossin(theta))*l
    
    def v(self, t):
        p = self.p
        e = self.e
        a = self.a
        u = self.u
        l = p/(1+e*ct)
        vel = sqrt(u*(2./l-1./a))
        return vel*self.prograde(t)
        
    @classmethod
    def from_config(cls, conf_parser, section, u):
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
        
