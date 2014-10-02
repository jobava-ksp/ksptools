from __future__ import division

import numpy as np
import numpy.linalg as npla
import scipy as sp
import scipy.constants as spconst
from scipy.integrate import odeint
from scipy.optimize import brentq, minimize, newton

def accel_func(T):
    def func(r, v, m):
        accel_g = -spconst.g
        accel_t = T/m
        accel_d = -8e-4*1.22*(v**2)*np.e**(-r/5e+3)
        return accel_g + accel_t + accel_d
    return func

def accel_func2d(T, fa):
    def func(r, v, m):
        rx, ry = r
        vx, vy = v
        ux = rx/npla.norm(r)
        uy = ry/npla.norm(r)
        accel_gx = -ux*spconst.g
        accel_gy = -uy*spconst.g
        accel_tx = (uy*np.cos(fa) - ux*np.sin(fa))*T/m
        accel_ty = (ux*np.cos(fa) + uy*np.sin(fa))*T/m
        
        vair_x = vx + uy*174
        vair_y = vy - ux*174
        vair = np.array([vair_x, vair_y])
        
        accel_dy = -8e-4*1.22*vair_y*(npla.norm(vair))*np.e**(-(npla.norm(r)-6e+5)/5e+3)
        accel_dx = -8e-4*1.22*vair_x*(npla.norm(vair))*np.e**(-(npla.norm(r)-6e+5)/5e+3)
        ay = accel_gy + accel_ty + accel_dy
        ax = accel_gx + accel_tx + accel_dx
        return np.array([ax, ay])
    return func

def dm_func(T, isp_0, isp_1):
    def func(r, v, m):
        isp = isp_0 + min(1, np.e**(-(npla.norm(r)-6e+5)/5e+3))*(isp_1 - isp_0)
        return -T/(isp*spconst.g)
    return func


def _burn_for_t(y0, t, af, df):
    afunc = lambda y, t: np.array(list(y[2:4]) + list(af(y[0:2], y[2:4], y[4])) + [df(y[0:2], y[2:4], y[4])])
    rx, ry, vx, vy, m1 = odeint(afunc, y0, np.arange(0, t, 0.01))[-1]
    return t, array([rx,ry]), array([vx, vy]), m1

def _burn_to_objective(y0, t0, af, df, objfunc):
    afunc = lambda y, t: np.array(list(y[2:4]) + list(af(y[0:2], y[2:4], y[4])) + [df(y[0:2], y[2:4], y[4])])
    def sim_func(t):
        rx, ry, vx, vy, m1 = odeint(afunc, y0, np.arange(0, t, 0.01))[-1]
        r = np.array([rx,ry])
        v = np.array([vx,vy])
        return objfunc(r, v, m1, t)
    t = brentq(sim_func, 0.5*t0, 2.5*t0)
    #t = minimize(lambda t: sim_func(t)**2, t0, method='CG').x[0]
    rx, ry, vx, vy, m1 = odeint(afunc, y0, np.arange(0, t, 0.01))[-1]
    return t, np.array([rx, ry]), np.array([vx, vy]), m1

def _initial_conditions_twr(r, v, Twr, T, isp_0, isp_1, fa):
    af = accel_func2d(T, fa)
    df = dm_func(T, isp_0, isp_1)
    m0 = T/(Twr*spconst.g)
    y0 = np.concatenate((list(r), list(v), [m0]))
    return af, df, m0, y0

def dv_to_alt(alt, Twr, T, isp_0, isp_1, r=[0,6e+5], v=[-174.0,0], fa=np.pi/2):
    af, df, m0, y0 = _initial_conditions_twr(r, v, Twr, T, isp_0, isp_1, fa)
    tf = np.sqrt((2*alt)/(spconst.g*(Twr-1)))
    def altfunc(r, v, m1, t):
        return (alt+6e+5) - npla.norm(r)
    t, r, v, m1 = _burn_to_objective(y0, tf, af, df, altfunc)
    idv = isp_0 * spconst.g * np.log(m0/m1)
    return t, r, v, m1, m0, idv

t, r, v, m1, m0, idv = dv_to_alt(9e+3, 2.7, 2130e+3, 318.5, 286.8)
print(t)
print('{} {}'.format(npla.norm(r), r))
print('{} {}'.format(npla.norm(v), v))
print((m1,m0,m0-m1))
print(idv)


