from __future__ import print_function, division

import collections
import scipy.optimize

from numpy import sqrt, dot, cos, sin, arccos, arcsin, pi, arccosh, sinh, tan, spacing, arctan, arctan2, finfo
from numpy.linalg import norm
from .util import arcvec, cossin

TransferParameters = collections.namedtuple("TransferParameters", ['u','r1','v1','r2','v2','t0','t1', 'dta'])

def commonplane(kep1, kep2, t0, t1):
    r1, v1 = kep1.rv(t0)
    r2, v2 = kep2.rv(t1)
    r2 = kep1.orientation / (kep2.orientation * r2)
    v2 = kep1.orientation / (kep2.orientation * v2)
    return r1, v1, r2, v2, lon_asc

def solve_transfer(u, kep1, kep2, t0, t1, method='ballistic', solver=SemiLatisSolver, output='kepler'):
    r1, v1 = kep1.rv(t0)
    r2, v2 = kep2.rv(t1)
    
    if method == 'ballistic':
        dta = arcvec(r1, r2)
        mnsol, err, solv1, solv2 = solve_transfer_planar(TransferParameters(u, r1, v1, r2, v2, t0, t1, dta), solver())
    
    if output == 'debug':
        return mnsol, err, solv1, solv2

def solve_transfer_planar(params, solver):
    x = solver.start(params)
    t0, t1 = params.t0, params.t1
    
    f = lambda i: (solver.tof(i, params) - (t1 - t0))**2
    fp = lambda i: 2*(solver.tof(i, params) - (t1 - t0))*solver.tofprime(i, params)
    
    p = scipy.optimize.minimize(f, x, jac=fp, bounds=[(solver.minvalue, solver.maxvalue)], method='TNC')
    err = sqrt(p.func)
    solv1, solv2 = solver.solution(p.x[0], params)
    return p, err, solv1, solv2


class TransferSolver(object):
    pass


class SemiLatisSolver(TransferSolver):
    def checkp(self, p, params):
        if p <= self.minvalue:
            print("bad p, {} < {}".format(p, self.minvalue))
        elif p >= self.maxvalue:
            print("bad p, {} > {}".format(p, self.maxvalue))
        #else:
        #    print("good p")
        
        #assert(p > self.minvalue)
        #assert(p < self.maxvalue)
        return p
    
    def checka(self, a, p, params):
        if abs(a) == float('inf'):
            print(a)
            return (finfo(a).max)**(1/3.)
        assert(a != 0)
        return a
    
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
        p = self.checkp(norm(p), params)
        
        r1, r2 = norm(params.r1), norm(params.r2)
        cosdt, sindt = self.cosdt, self.sindt
        
        f = 1 - r2*(1-cosdt)/p
        g = r1*r2*sindt/sqrt(params.u*p)
        fp = sqrt(params.u/p)*tan(self.theta/2)*((1-cosdt)/p - 1/r1 - 1/r2)
        gp = 1 - r1*(1-cosdt)/p
        
        m, k, l = self.m, self.k, self.l
        a = m*k*p/((2*m-l**2)*p**2 + 2*k*l*p - k**2)
        
        a = self.checka(a,p,params)
        
        if a > 0:
            cosE = 1 - r1*(1-f)/a
            sinE = -r1*r2*fp/sqrt(params.u*a)
            #E = arctan(sinE/cosE)
            #E = arcsin(cosE)
            E = arccos(cosE)
            if sinE < 0:
               E = 2*pi - E
            assert(E >= 0)
            return g + sqrt(a**3/params.u)*(E-sinE)
        elif a < 0:
            coshF = 1 - r1*(1-f)/a
            F = arccosh(coshF)
            return g + sqrt((-a)**3/params.u)*(sinh(F)-F)
        else:
            print(a)
            raise NotImplementedError 
        
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
        #n = 3
        #eps = spacing(n*p)
        return scipy.optimize.approx_fprime(p, lambda i: self.tof(i, params), 1e-4)
        
