
from numpy import arccos, arcsin, arctan, arctan2, array, cross, cos, dot, sin, sqrt
from numpy.linalg import norm
from scipy.optimize import newton

def project(a,b):
    return (dot(a,b)/dot(b,b))*b

def reject(a,b):
    return a - (dot(a,b)/dot(b,b))*b


class KeplerOrbit(object):
    def __init__(self, body, p, e, inc, lon_asc, arg_pe, M, epoch):
        self.body = body
        self.p = p
        self.e = e
        self.a = p/(1-e**2)
        self.inc = inc
        self.lon_asc = lon_asc
        self.arg_pe = arg_pe
        self.epoch = epoch
        #self.true_anom = theta
        #self.E_anom = arccos((e + cos(theta))/(1 + e*cos(theta)))
        self.mean_motion = sqrt(body.std_g_param/self.a**3)
        #self.mean_anom = self.E_anom - e*sin(self.E_anom)
        self.mean_anom = M
     
    @classmethod
    def from_planet_paremters(cls, body, a, e, i, arg_pe, lon_asc, M, epoch):
        p = (1-e**2)*a
        return cls(body, p, e, i, lon_asc, arg_pe, M, epoch)
    
    @classmethod
    def from_rvbody(cls, r, v, body, epoch):
        p,e,t,inc,lon_asc,arg_pe = KeplerOrbit._from_rvu(r,v,body.std_g_param)
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
        x = cos(theta)*(r/norm(r)) - sin(theta)*(vt/norm(vt))
        y = sin(theta)*(r/norm(r)) + cos(theta)*(vt/norm(vt))
        z = cross(x,y)

        E = arccos((e+cos(theta))/(1+e*cos(theta)))
        M = E - e*sin(E)
        
        lon_asc = arctan2(z[0], -z[1])
        inc = arctan2(sqrt(z[0]**2+z[1]**2), z[2])
        arg_pe = arctan2(x[2],y[2])

        return p,e,M,inc,lon_asc,arg_pe

    def to_rvbody(self, t):
        dt = t - self.epoch
        M = self.mean_anom + (self.mean_motion * t)
        f = lambda nE: nE - e*np.sin(nE) - M - n*t
        E = newton(f,0)
        u = self.body.std_g_param
        a = self.a
        p = self.p
        e = self.e
        theta = arctan(sqrt(((1+e)*tan(E/2)**2)/(1-e)))/2
        len_r = (a*(1-e**2))/(1+e*cos(theta))
        O, i, w = self.lon_asc, self.inc, self.arg_pe
        x = array([
                cos(O)*cos(w) - sin(O)*cos(i)*sin(w),
                sin(O)*cos(w) + cos(O)*cos(i)*sin(w),
                sin(i)*sin(w)])
        y = array([
                -cos(O)*sin(w) - sin(O)*cos(i)*cos(w),
                -sin(O)*sin(w) + cos(O)*cos(i)*cos(w),
                sin(i)*cos(w)])
        z = array([
                sin(i)*sin(O),
                -sin(i)*cos(O),
                cos(i)])
        unit_r = cos(theta)*x + sin(theta)*y
        unit_t = cos(theta)*y - sin(theta)*x
        r = len_r*unit_r
        vel = sqrt(u*(2/len_r-1/a))
        vt = (sqrt(u*p)/len_r)*unit_t
        vr = (vel-norm(vt))*unit_r
        return r, vt + vr, body

