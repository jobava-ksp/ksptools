from numpy import array, concatenate, dot, zeros
from numpy.linalg import norm
from scipy.integrate import odeint

def simrv(stv0, t0, t1, accel=lambda r,v,t: zeros(3)):
    rv0 = concatenate((stv0.r, stv0.v))
    def dfunc(y, t):
        return concatenate((y[3:6], accel(y[0:3], y[3:6] ,t)))
    rv_array = odeint(dfunc, rv0, [t0, t1])
    return type(stv0)(rv_array[-1][0:3], rv_array[-1][3:6])

def gravity_accel_func(u):
    def gafunc(r, v, t):
        return -(r/norm(r))*(u/dot(r,r))
    return gafunc


