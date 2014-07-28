from numpy import arccos, arcsin, arctan, arctan2, array, cross, cos, dot, mat, sin, sqrt, tan
from numpy.linalg import norm
from scipy.optimize import newton

from . import unit

def project(a,b):
    return (dot(a,b)/dot(b,b))*b

def reject(a,b):
    return a - (dot(a,b)/dot(b,b))*b

def rotx(t):
    return mat([[1,       0,      0],
                [0,  cos(t), sin(t)],
                [0, -sin(t), cos(t)]])

def roty(t):
    return mat([[cos(t), 0, -sin(t)],
                [0,      1,       0],
                [sin(t), 0,  cos(t)]])

def rotz(t):
    return mat([[cos(t), -sin(t), 0],
                [sin(t),  cos(t), 0],
                [     0,       0, 1]])
                

class EulerAngle(object):    
    def __init__(self, p, t, s):
        Ap = rotz(p)  
        At = rotx(t)
        As = rotz(s)
        
        self.phi, self.theta, self.psi = p, t, s
        self.Ap, self.At, self.As = Ap, At, As
        self.A = Ap*At*As
    
    @classmethod
    def from_pts(cls,p,t,s):
        return cls(p,t,s)
    
    def rotate(self, r):
        return (self.A*(mat(r).T)).A1

class KeplerOrbit(object):
    def __init__(self, body, a, e, inc, lon_asc, arg_pe, M, epoch):
        self.body = body
        self.e = e
        self.a = a
        self.orient = EulerAngle.from_pts(lon_asc, inc, arg_pe)
        #self.inc = inc
        #self.lon_asc = lon_asc
        #self.arg_pe = arg_pe
        self.epoch = epoch
        self.mean_motion = sqrt(body.std_g_param/self.a**3)
        self.mean_anom = M
     
    @classmethod
    def from_planet_paremters(cls, body, a, e, i, arg_pe, lon_asc, M, epoch):
        #p = (1-e**2)*a
        return cls(body, a, e, i, lon_asc, arg_pe, M, epoch)
    
    @classmethod
    def from_rvbody(cls, r, v, body, epoch, u2=0):
        p,e,t,inc,lon_asc,arg_pe = KeplerOrbit._from_rvu(r,v,body.std_g_param + u2)
        return cls(body, p, e, inc, lon_asc, arg_pe, t, epoch)

    @staticmethod
    def _from_rvu(r, v, u):
        vr = project(v,r)
        vt = reject(v,r)
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

        return p,e,M,inc,lon_asc,arg_pe

    #def rv(self, t, u2=0):
    #    u = self.body.std_g_param + u2
    #    a = self.a
    #    p = self.p
    #    e = self.e
    #    
    #    theta = self.true_anom(t)
    #    #x,y,z = self.xyz()
    #    
    #    len_r = p/(1+e*cos(theta))
    #    #len_r = (a*(1-e**2))/(1+e*cos(theta))
    #    unit_r = cos(theta)*x + sin(theta)*y
    #    unit_t = cos(theta)*y - sin(theta)*x
    #    r = len_r*unit_r
    #    vel = sqrt(u*(2/len_r-1/a))
    #    vt = sqrt(u*p)/len_r
    #    vr = sqrt(vel**2-vt**2)
    #    return r, vt*unit_t + vr*unit_r

    def rv(self, t, u2=0):
        u = self.body.std_g_param + u2
        a = self.a
        e = self.e
        b = a*sqrt(1-e**2)
        p = a*(1-e**2)
        ### solve in the orbital plane ###
        theta = self.true_anom(t)
        r = p*array([cos(theta), sin(theta), 0])/(1+e*cos(theta))
        v = array([-r[1]/b**2, r[0]/a**2, 0])
        ### rotate into equitorial plane ###
        r = self.orient.rotate(r)
        v = self.orient.rotate(v)
        return r,v
    
    def E_anom(self, t):
        dt = t - self.epoch
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
    
    #def theta(self, t):
    #    return self.true_anom(t) + self.orient.phi + self.orient.psi
    
    #def xyz(self):
    #    O, i, w = self.lon_asc, self.inc, self.arg_pe
    #    x = array([
    #            cos(O)*cos(w) - sin(O)*cos(i)*sin(w),
    #            sin(O)*cos(w) + cos(O)*cos(i)*sin(w),
    #            sin(i)*sin(w)])
    #    y = array([
    #            -cos(O)*sin(w) - sin(O)*cos(i)*cos(w),
    #            -sin(O)*sin(w) + cos(O)*cos(i)*cos(w),
    #            sin(i)*cos(w)])
    #    z = array([
    #            sin(i)*sin(O),
    #            -sin(i)*cos(O),
    #            cos(i)])
    #    return x,y,z
        
    @classmethod
    def from_config(cls, conf_parser, section, body):
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
        return cls.from_planet_paremters(body, a, e, i, arg_pe, lon_asc, M, epoch)
        
