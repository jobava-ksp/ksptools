from .._vector import statevector

from numpy import array, concatenate, dot, zeros
from numpy.linalg import norm
from scipy.integrate import odeint

import scipy.constants as spconst


def simrv(stv0, t0, t1, accel=lambda r, v, t: zeros(3)):
    rv0 = concatenate((stv0.r, stv0.v))

    def dfunc(y, t):
        r = y[0:3]
        v = y[3:6]
        a = accel(r, v, t)
        return concatenate((v, a))

    rv_array = odeint(dfunc, rv0, [t0, t1])
    return type(stv0)(rv_array[-1][0:3], rv_array[-1][3:6])


def simrvm(stv0, m0, t0, t1, accel=lambda r, v, t, m: zeros(3), dmfunc=lambda r, v, t, m: 0):
    rv0 = concatenate((stv0.r, stv0.v, array([m0])))

    def dfunc(y, t):
        r = y[0:3]
        v = y[3:6]
        m = y[6]
        a = accel(r, v, t, m)
        dm = dmfunc(r, v, t, m)
        return concatenate((v, a, array([dm])))

    rv_array = odeint(dfunc, rv0, [t0, t1])
    return statevector(rv_array[-1][0:3], rv_array[-1][3:6]), rv_array[6]


def gravity_accel_func(u):

    def gafunc(r, v, t, *args):
        return -(r/norm(r))*(u/dot(r,r))
    return gafunc


def thrust_accel_func(d, isp, dm):

    def tafunc(r, v, t, m, *args):
        return d*(spconst.g*isp*dm)/m
    
    def dmfunc(r, v, t, m, *args):
        return -dm
    
    return tafunc, dmfunc

