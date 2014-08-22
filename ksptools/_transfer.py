from __future__ import print_function, division

import collections
import scipy.optimize

from numpy import array, dot, cross, sqrt, pi
from numpy import cos, sin, tan, arccos, arcsin, arctan, arctan2
from numpy import sinh, arccosh
from numpy import spacing, finfo
from numpy.linalg import norm, solve

from .util import arcvec, cossin, unit, rotz, rotx, rotvec, Ax
from .orbit import KeplerOrbit, Patch

TransferParameters = collections.namedtuple("TransferParameters", ['u','r1','v1','r2','v2','t1','t2', 'dta'])


class TransferSolver(object):
    def start(self, params):
        raise NotImplementedError
    
    def tof(self, x, params):
        raise NotImplementedError
    
    def tofprime(self, x, params):
        raise NotImplementedError
    
    def solution(self, x, params):
        raise NotImplementedError


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
        if output == 'pass':
            return mnsol, err, solv1, solv2
        elif output == 'kepler':
            return [(KeplerOrbit.from_rvu(r1, solv1, u, t1), t1, t2)]
    elif method == 'planechange':
        R, Rinv, delta_ta1, inc = solve_planechange(u, r1, v1, r2, v2)
        mnsol, err, solv1, solv2 = solve_transfer(u, r1, v1, Ax(R, r2), Ax(R, v2), t1, t2, method='ballistic', solver=solver, output='pass')
        kep1 = KeplerOrbit.from_rvu(r1, solv1, u, t1)
        kep2 = KeplerOrbit.from_rvu(r2, Ax(Rinv, solv2), u, t2)
        tc = kep1.time_at_true_anomaly(delta_ta1 + kep1.true_anomaly(t1))
        if output == 'pass':
            return mnsol, err, solv1, solv2, kep2.velocity(tc) - kep1.velocity(tc)
        if output == 'kepler':
            return [(kep1, t1, tc), (kep2, tc, t2)]

def solve_transfer_planar(params, solver):
    x = solver.start(params)
    t1, t2 = params.t1, params.t2
    
    f = lambda i: (solver.tof(i, params) - (t2 - t1))**2
    fp = lambda i: 2*(solver.tof(i, params) - (t2 - t1))*solver.tofprime(i, params)
    
    p = scipy.optimize.minimize(f, x, jac=fp, bounds=[(solver.minvalue, solver.maxvalue)], method='TNC')
    err = sqrt(p.fun)
    solv1, solv2 = solver.solution(p.x[0], params)
    return p, err, solv1, solv2

def solve_planechange(u, r1, v1, r2, v2):
    normal1 = unit(cross(r1,v1))
    normal2 = unit(cross(r2,v2))
    axis = cross(normal1,normal2)
    ta = arcvec(r1, axis)
    sint = norm(cross(normal1,normal2)/(norm(normal1)*norm(normal2)))
    R = rotvec(axis, -arcsin(sint))
    Rinv = rotvec(axis, arcsin(sint))
    return R, Rinv, ta, arcsin(sint)


def solve_flyby(u, distance, vsoi, soi, eta, insertion=False):
    ## Velocity at periapsis ##
    vp = sqrt(norm(vsoi)**2 + 2*u*((norm(soi)-distance)/(norm(soi)*distance)))
    
    ## Basic Orbital Parameters ##
    a = 1/(2/distance-(vp**2)/u)     # semi major axis
    e = (distance*vp**2)/u - 1       # eccentricity
    cost = (a*(1-e**2)-soi)/(e*soi)  # cosine of true anomaly passing soi
    t = arccos(cost)                 # true anomaly (assuming an ejection)
    if insertion:                    # fix for insertion
        t = -t
    sint = sin(t)                    # sine of true anomaly
    
    ## Orientation Parameters ##
    fa = arctan(e*sint/(1+e*cost))
    fpath = 0.5*pi + t - fa
    
    vi, vj, vk = vsoi                                  # velocity by component
    vx, vy = norm(vsoi) * array([cos(fpath), sin(fpath)])
    
    if vk >= 0: # ascending
        inc = arcsin(vk/vy)              # inclenation for ascending node at pe
        argpe = 0                        # argument pe (right on top of ascending node)
        cos_lonasc, sin_lonasc = solve(array([[vx, -cos(inc)*vy],[cos(inc)*vy, vx]]), array([vi, vj]))
        lonasc = arctan2(sin_lonasc, cos_lonasc)
    else: # descending
        inc = arcsin(-vk/vy)             # inclenation for descending node at pe
        argpe = pi                       # argument pe (opposite of ascending node)
        cos_lonasc, sin_lonasc = solve(array([[-vx, cos(inc)*vy],[-cos(inc)*vy, -vx]]), array([vi, vj]))
        lonasc = arctan2(sin_lonasc, cos_lonasc)
    A = rotz(lonasc)*rotx(inc)*rotz(argpe)       # rotation matrix
    
    ## Position at Soi ##
    rplanar = array([cost, sint, 0])*soi  # planar position
    rsoi = Ax(A, rplanar)                 # position at soi
    
    ## Position at Pe
    rplanar = array([1, 0, 0])*distance   # planar position
    rpe = Ax(A, rplanar)                  # position at pe
    
    ## Time offset ##
    coshF = (e + cost)/(1+e*cost)   # cosine hyperbolic of eccentric anomaly
    F = arccosh(coshF)              # eccentric anomaly (assuming ejection)
    mm = sqrt(u/(-a)**3)            # mean motion
    delta_time = (e*sinh(F)-F)/mm   # delta time

    if insertion:
        return vp, rpe, rsoi, eta - delta_time
    else:
        return vp, rpe, rsoi, eta + delta_time

