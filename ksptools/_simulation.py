from scipy.integrate import romberg
from numpy import dot
from numpy.linalg import norm

from .util import unit, Ax

class Controller(object):
    pass


def accel_g(u, r):
    return -(u*unit(r))/dot(r,r)

def accel_d(cf, pfunc, r, v):
    return -unit(v)*0.5*8e-3*cf*pfunc(r)*dot(v,v)

def accel_f(D, engines, atmfunc, mass0, t0, t, r, v):
    dt = t - t0
    thrust = sum(e.thrust for e in engines)
    ff = sum(e.thrust/e.isp(atmfunc(r)) for e in engines)
    return Ax(D, unit(r)) * (thrust/(mass0 - dt*ff))

def runpass(r0, v0, t0, mass

def run(initial_state, controller):
    r0 = state.position
    v0 = state.velocity
    t0 = state.epoch
    
