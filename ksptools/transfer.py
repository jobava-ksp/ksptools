from __future__ import print_function, division

import collections
import scipy.optimize

from numpy import sqrt, dot, cos, sin, arccos, arcsin, pi, arccosh, sinh, tan, spacing, arctan, arctan2, finfo
from numpy.linalg import norm

from .util import arcvec, cossin
from .orbit import KeplerOrbit

TransferParameters = collections.namedtuple("TransferParameters", ['u','r1','v1','r2','v2','t1','t2', 'dta'])


class TransferSolver(object):
    pass


class SemiLatisSolver(TransferSolver):
    def checkpa(self, p, params):
        m, k, l = self.m, self.k, self.l
        
        under = ((2*m-l**2)*p**2 + 2*k*l*p - k**2)
        over = m*k*p
        
        if abs(under) > 0:
            a = over/under
        elif over == 0:
            a = 0
        else:
            a = (finfo(over/(under+spacing(over))).max)**(1/3.)
        
        if p < spacing(params.u):
            p = spacing(params.u)
        
        return p, a
    
    #def checka(self, a, p, params):
    #    if abs(a) == float('inf'):
    #        #print(a)
    #        return (finfo(a).max)**(1/3.)
    #    #assert(a != 0)
    #    if a is None:
    #        return 0
    #    return a
    
    def start(self, params):
        cosdt, sindt, _ = cossin(params.dta)
        self.cosdt = cosdt
        self.sindt = sindt
        self.theta = params.dta
        
        r1, r2 = norm(params.r1), norm(params.r2)
        self.k = r1*r2*(1-cosdt)
        self.m = r1*r2*(1+cosdt)
        self.l = r1 + r2
        
        p_i = self.k/(self.l + sqrt(2*self.m))
        p_ii = self.k/(self.l - sqrt(2*self.m))
        
        if self.theta < pi:
            self.minvalue = p_i
            self.maxvalue = float('inf')
            return p_i*3
        elif self.theta > pi:
            self.minvalue = 0
            self.maxvalue = p_ii
            return p_ii/2
        else:
            raise Exception
    
    def tof(self, p, params):
        p, a = self.checkpa(norm(p), params)
        
        r1, r2 = norm(params.r1), norm(params.r2)
        cosdt, sindt = self.cosdt, self.sindt
        
        f = 1 - r2*(1-cosdt)/p
        g = r1*r2*sindt/sqrt(params.u*p)
        fp = sqrt(params.u/p)*tan(self.theta/2)*((1-cosdt)/p - 1/r1 - 1/r2)
        gp = 1 - r1*(1-cosdt)/p
        
        if a > 0:
            cosE = 1 - r1*(1-f)/a
            sinE = -r1*r2*fp/sqrt(params.u*a)
            #E = arctan(sinE/cosE)
            #E = arcsin(cosE)
            _a = min(1.0, max(-1.0, cosE))
            E = arccos(_a)
            if sinE < 0:
               E = 2*pi - E
            assert(E >= 0)
            return g + sqrt(a**3/params.u)*(E-sinE)
        elif a < 0:
            coshF = 1 - r1*(1-f)/a
            F = arccosh(coshF)
            if F < 0:
                F = 2*pi + F
            return g + sqrt((-a)**3/params.u)*(sinh(F)-F)
        elif a == 0:
            return 0
        
    def solution(self, p, params):
        r1, r2 = norm(params.r1), norm(params.r2)
        rv1, rv2 = params.r1, params.r2
        cosdt, sindt = self.cosdt, self.sindt
        
        f = 1 - r2*(1-cosdt)/p
        g = r1*r2*sindt/sqrt(params.u*p)
        fp = sqrt(params.u/p)*tan(self.theta/2)*((1-cosdt)/p - 1/r1 - 1/r2)
        gp = 1 - r2*(1-cosdt)/p
        
        v1 = (rv2 - f*rv1)/g
        v2 = fp*rv1 + gp*v1
        
        return v1, v2
    
    def tofprime(self, p, params):
        eps = spacing(params.u)
        #eps = 5.0e-5
        return scipy.optimize.approx_fprime(p, lambda i: self.tof(i, params), eps)


def solve_transfer(u, r1, v1, r2, v2, t1, t2, method='ballistic', solver=SemiLatisSolver, output='kepler'):
    if method == 'ballistic':
        dta = arcvec(r1, r2)
        mnsol, err, solv1, solv2 = solve_transfer_planar(TransferParameters(u, r1, v1, r2, v2, t1, t2, dta), solver())
        if output == 'debug':
            return mnsol, err, solv1, solv2
        elif output == 'kepler':
            return KeplerOrbit.from_rvu(r1, solv1, u, t1)

def solve_transfer_planar(params, solver):
    x = solver.start(params)
    t1, t2 = params.t1, params.t2
    
    f = lambda i: (solver.tof(i, params) - (t2 - t1))**2
    fp = lambda i: 2*(solver.tof(i, params) - (t2 - t1))*solver.tofprime(i, params)
    
    p = scipy.optimize.minimize(f, x, jac=fp, bounds=[(solver.minvalue, solver.maxvalue)], method='TNC')
    err = sqrt(p.fun)
    solv1, solv2 = solver.solution(p.x[0], params)
    return p, err, solv1, solv2


def solve_flyby(u, rp, vel0, soi):
    vp = sqrt(norm(v0)**2 + 2*u/rp) # orbital speed at pe
    a = 1/(2/rp-(vp**2)/u)          # semi major axis
    e = (rp*vp**2)/u - 1            # eccentricity
    cost = (a*(1-e**2)-soi)/(e*soi) # cosine of true anomaly passing soi
    vel0k = project(vel0, unitk)[2] # velocity going up
    
    coshF = (e + cost)/(1+e*cost)   # cosine hyperbolic of eccentric anomaly
    F = arccosh(coshF)              # eccentric anomaly (without sign)
    mm = sqrt((-a)**3/u)            # mean motion
    delta_time = F/mm               # delta time
    return vp, a, e, cost, F, delta_time, vel0k

def solve_ejection(u, rp, vel0, soi):
    vp, a, e, cost, F, delta_time, vel0k = solve_flyby(u, rp, vel0, soi)
    
    dt = a*(e**2-1)*cost/(1+e*cost) # dt for dx/dt and dy/dt (planar velocity)
    vplanar_y = norm(vel0)*dtdy     # y component of planar velocity at soi
    
    if vel0k >= 0: # ascending
        inc = arcsin(vel0k/vplanar_y) # inclenation for ascending node at pe
    else: # descending
        inc = arcsin(-vel0k/vplanar_y) # inclenation for descending node at pe
    
