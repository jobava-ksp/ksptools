from __future__ import print_function, division

import collections
import scipy.optimize

from numpy import sqrt

from .orbit import spec_energy, time_from_pe, semi_major_axis, eccentricity
from .util import arcvec

TransferParamters = collections.namedtuple("TransferParamters", ['u','r1','v1','r2','v2','t0','t1', 'dta'])

def commonplane(kep1, kep2, t0, t1):
    r1, v1 = kep1.rv(t0)
    r2, v2 = kep2.rv(t1)
    r2 = kep1.orientation / (kep2.orientation * r2)
    v2 = kep1.orientation / (kep2.orientation * v2)
    return r1, v1, r2, v2, lon_asc

def solve_transfer(u, kep1, kep2, t0, t1, method='ballistic'):
    r1, v1 = kep1.rv(t0)
    r2, v2 = kep2.rv(t1)
    
    if method == 'ballistic':
        dta = arcvec(r1, r2)
        return slove_transfer_planar(TransferParameters(u, r1, v1, r2, v2, t0, t1, dta))

def solve_transfer_planar(params, pick, toffunc, compute):
    x = pickv(params)
    scipy.optimize.minimize(lambda i: abs(toffunc(i, params) - (t1 - t0)), x)
    return compute(x)

def pickv(params):
    from numpy.linalg import norm
    return params.v1

def tofv(v0, params):
    u, r1 = params.u, params.r1
    a = semi_major_axis(u=u, r=r1, v=v0)
    e = eccentricity(u=u, a=a, r=r1, v=v0)
    t1 = time_from_pe(u=u, a=a, e=e, r=r1)
    t2 = time_from_pe(u=u, a=a, e=e, r=r2)
    return t2 - t1
