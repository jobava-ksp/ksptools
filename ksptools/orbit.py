from __future__ import division

class KeplerOrbit(object):
    def __init__(self, u, a, e, inc, lon_asc, arg_pe, M, epoch):
        from numpy import sqrt, pi
        from .util import EulerAngle
        self.u = u
        self.e = e
        self.p = a*(1-e**2)
        self.a = a
        self.orient = EulerAngle.from_pts(lon_asc, inc, arg_pe)
        self.mean_motion = sqrt(u/self.a**3)
        #while M < 0:
        #    M += 2*pi
        self.mean_anom = M
     
    @classmethod
    def from_planet_paremters(cls, u, a, e, i, arg_pe, lon_asc, M, epoch):
        return cls(u, a, e, i, lon_asc, arg_pe, M, epoch)

    #@classmethod
    #def from_rvu(cls, r, v, u, epoch=0.):
        #from numpy import arccos, arctan, arctan2, cross, cos, sin, sqrt
        #from numpy.linalg import norm
        #from util import project, reject
        #vr = project(v,r)
        #vt = reject(v,r)
        #p = (norm(r)*norm(vt))**2/u
        #e = sqrt(1 + 2*(norm(v)**2/2 - u/norm(r))*(norm(cross(r,v)))**2/u**2)
        #if e == 0.0:
        #    theta = 0
        #else:
        #    v0 = sqrt(u/p)
        #    acos_arg = min(1,max(-1,norm(vt)/(v0*e)-1/e))
        #    theta = arccos(acos_arg)
        #E = arccos((e+cos(theta))/(1+e*cos(theta)))
        #M = E - e*sin(E)
        #
        #x = cos(theta)*(r/norm(r)) - sin(theta)*(vt/norm(vt))
        #y = sin(theta)*(r/norm(r)) + cos(theta)*(vt/norm(vt))
        #z = cross(x,y)
        #lon_asc = arctan2(z[0], -z[1])
        #inc = arctan2(sqrt(z[0]**2+z[1]**2), z[2])
        #arg_pe = arctan2(x[2],y[2])

        ##return p,e,M,inc,lon_asc,arg_pe
        #return cls(u, p, e, inc, lon_asc, arg_pe, M, epoch)
    
    @classmethod
    def from_rvu(cls, rvec, velvec, u, epoch=0.):
        from numpy import cross, array, dot, arccos, sin, cos, pi
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
        
        argpe = arccos(dot(n,evec)/(norm(n)*norm(evec)))
        if evec[2] > 0: argpe = 2*pi - argpe
        
        costa = dot(evec,rvec)/(norm(evec)*norm(rvec))
        E = arccos((e+costa)/(1+e*costa))
        M = E - e*sin(E)
        if dot(rvec, velvec) > 0: M = 2*pi - M
        
        return cls(u, a, e, i, la, argpe, M, epoch)

    def E_anom(self, dt):
        from numpy import cos, sin
        from scipy.optimize import newton
        M = self.mean_anom + (self.mean_motion * dt)
        e = self.e
        f = lambda nE: nE - e*sin(nE) - M
        fp = lambda nE: 1 - e*cos(nE)
        fp2 = lambda nE: e*sin(nE)
        return newton(f,0,fp,fprime2=fp2)
    
    def true_anom(self, t):
        from numpy import arctan, tan, sqrt
        E = self.E_anom(t)
        e = self.e
        return 2*arctan(sqrt((1.+e)/(1.-e))*tan(E/2.))
    
    def flight_angle(self, t, theta=None):
        from .util import cossin
        from numpy import arctan
        e = self.e
        if theta is None:
            theta = self.true_anom(t)
        ect, est, _ = e*cossin(theta)
        return -arctan(est/(1.+ect))
    
    def prograde(self, t, theta=None):
        from .util import cossin, Ax, rotz
        from numpy import array
        e = self.e
        if theta is None:
            theta = self.true_anom(t)
        ct, st, _ = cossin(theta)
        
        pgd = Ax(rotz(self.flight_angle(t, theta)+theta), array([0., 1., 0.]))
        #pgd = Ax(rotz(theta), array([0., 1., 0.]))
        return self.orient.rotate(pgd)
    
    def radialin(self, t, theta=None):
        from .util import cossin, Ax, rotz
        from numpy import array
        e = self.e
        if theta is None:
            theta = self.true_anom(t)
        ct, st, _ = cossin(theta)
        rad = Ax(rotz(self.flight_angle(t, theta)), array([-1., 0., 0.]))
        return self.orient.rotate(rad)
    
    def normal(self, t):
        from numpy import array
        return self.orient.rotate(array([0.,0.,1.]))
    
    def r(self, t, theta=None):
        from .util import cossin
        from numpy import cos
        p,e = self.p, self.e
        if theta is None:
            theta = self.true_anom(t)
        l = p/(1+e*cos(theta))
        return self.orient.rotate(cossin(theta))*l
    
    def v(self, t, theta=None):
        from .util import cossin, Ax, rotz
        from numpy import array, arctan2, sqrt
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
        from numpy import pi
        return 2*pi/self.mean_motion
        
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
        
