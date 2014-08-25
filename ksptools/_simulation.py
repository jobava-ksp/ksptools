from scipy.integrate import romberg
from numpy import dot
from numpy.linalg import norm

from .util import unit, Ax

class Controller(object):
    pass


def getatmfunc(body):
    if body.atmposhere is None:
        return lambda x: 0
    else:
        return lambda r: body.atmposhere.atm(norm(r)-body.eq_radius)

def getpfunc(body):
    if body.atmposhere is None:
        return lambda x: 0
    else:
        return lambda r: body.atmposhere.p(norm(r)-body.eq_radius)

def accel_g(u, r):
    return -(u*unit(r))/dot(r,r)

def accel_d(cf, pfunc, r, v):
    return -unit(v)*0.5*8e-3*cf*pfunc(r)*dot(v,v)

def accel_f(d, engines, atmfunc, mass, r, v):
    thrust = sum(e.thrust for e in engines)
    ff = sum(norm(e.thrust)/e.isp(atmfunc(r)) for e in engines)
    return d * (thrust/mass)

def physics_pass(

def run(initial_state, controller):
    r = state.position
    v = state.velocity
    t = state.epoch
    u = state.refbody.std_g_param
    atmfunc = getatmfunc(state.refbody)
    pfunc = getpfunc(state.refbody)
    
