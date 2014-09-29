from .._vector import statevector
from .._math import unit

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


def _thrust_accel_space_func(d, isp, dm):
    def tafunc(r, v, t, m, *args):
        return d(r,v)*(spconst.g*isp*dm)/m
    def dmfunc(r, v, t, m, *args):
        return -dm
    return tafunc, dmfunc


def _thrust_accel_atm_func(d, dm, atm, isp_1, isp_0):
    def tafunc(r, v, t, m, *args):
        _, a, _, _ = atm.atm_state(statevector(r,v),t)
        a = min(1, max(0, a))
        isp = isp_0 + a*(isp_1 - isp_0)
        return d(r,v)*(spconst.g*isp*dm)/m
    def dmfunc(r, v, t, m, *args):
        return -dm
    return tafunc, dmfunc


def thrust_accel_func(d, isp, dm, atm=None, isp_atm=None):
    if atm is not None:
        return _thrust_accel_atm_func(d, dm, atm, isp_atm, isp)
    else:
        return _thrust_accel_space_func(d, isp, dm)


def drag_acccel_func(atm, coefd=0.02):
    def dafunc(r, v, t, m, *args):
        p, _, _, v_air = atm.atm_state(statevector(r,v), t)


class Controller(object):
    def __init__(self):
        pass
    
    @classmethod
    def _accel_g(cls, r, body):
        return -(unit(r))*(body.GM/dot(r,r))
    
    @classmethod
    def _accel_thrust(cls, Tdird, T, m):
        return Tdir*T/m
    
    @classmethod
    def _accel_drag(cls, p, v_air, coefd):
        return -0.5*unit(v_air)*dot(v_air, v_air)*p*coefd*0.008
    
    @classmethod
    def _thrust_state(cls, a, engine_state):
        if engin_state is not None:
            Tmax, isp_0, isp_1, throttle = engine_state
            isp = isp_0 + min(1, max(0, a))*(isp_1 - isp_0)
            Tmax = isp_0 * ff * spconst.g
            T = Tmax * throttle
            ff = T/(isp * spconst.g)
            return T, ff
        else:
            return 0, 0
    
    def _update(self, m_rv, body_ijk_llav, pav, t):
        raise NotImplementedError
    
    def _do_state(self, r, v, m, body, t):
        cls = type(self)
        lat, lon, alt, v_air = body.llav(statevector(r,v), t)
        si, sj, sk = body.ijk_by_lla(lat, lon, alt)
        p, atm = body.atmstate_by_alt(alt, t)
        return (m, r, v), (body, si, sj, sk, lat, lon, alt), (p, atm, v_air), t
    
    def _do_control(self, m_rv, body_ijk_lla, pav, t):
        return self._update(m_rv, body_ijk_lla, pav, t)
    
    def _accel_dm(self, state, control):
        cls = type(self)
        ## state ##
        m_rv, body_ijk_llav, pav, t = state
        m, r, v = m_rv
        body, _, _, _, _, _, alt = body_ijk_lla
        p, a, v_air = pav
        
        ## control ##
        engine_state, coefd, Tdir = control
        T, ff = cls._thrust_state(a, engine_state)
        
        accel = cls._accel_g(r,body) + cls._accel_thrust(Tdir, T, m) + cls._accel_drag(p, v-v_air, coefd)
        return accel, ff
    
    def _ode_func(self, body):
        def func(self, y, t):
            r = y[0:3]
            v = y[3:6]
            m = y[6]
            state = self._do_state(r,v,m,body,t)
            control = self._do_control(*state)
            a, dm = self._accel_dm(state, control)
            return concatenate((v,a,[dm]))
        return func
    
    def _step(self, r, v, m, t, dt, odefunc):
        rvm0 = concatenate((r,v,[m]))
        rvm1 = odeint(odefunc, rvm0, [t, t+dt])[-1]
        return rvm1[0:3], rvm1[3:6], rvm1[6]
    
    def sim(self, stv, m, t, dt, body, pred=None):
        r, v = stv.r, stv.v
        func = self._ode_func(body)
        if pred is None:
            r, v, m = self._step(r, v, m, t, dt, func)
            t += dt
        else:
            while pred(r, v, m, t):
                r, v, m = self._step(r, v, m, t, dt, func)
                t += dt
        return statevector(r,v), m
    
    
