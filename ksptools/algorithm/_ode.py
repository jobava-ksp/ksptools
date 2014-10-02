from .._vector import statevector
from .._math import unit

from numpy import array, arange, concatenate, dot, zeros
from numpy.linalg import norm
from scipy.integrate import odeint
from scipy.optimize import brentq

import scipy.constants as spconst


class Controller(object):
    def __init__(self):
        pass
    
    def _initialize(self):
        self._rec = False
        self._record_log = None
    
    def _finalize(self):
        del self._rec
        del self._record_log
    
    @classmethod
    def _accel_g(cls, r, body):
        return -(unit(r))*(body.GM/dot(r,r))
    
    @classmethod
    def _accel_thrust(cls, Tdir, T, m):
        return (Tdir*T)/m
    
    @classmethod
    def _accel_drag(cls, p, v_air, coefd, m):
        return -(0.5)*unit(v_air)*p*dot(v_air, v_air)*coefd*8.0e-3
    
    @classmethod
    def _thrust_state(cls, a, engine_state):
        if engine_state is not None:
            Tmax, isp_0, isp_1, throttle = engine_state
            isp = isp_0 + min(1, max(0, a))*(isp_1 - isp_0)
            T = Tmax * throttle
            ff = T/(isp * spconst.g)
            return T, ff
        else:
            return 0, 0
    
    def _record_state(t, *args):
        self._record_log += [(t, args)]
    
    def _update(self, m_rv, body_ijk_llav, pav, t):
        raise NotImplementedError
    
    def _do_state(self, r, v, m, body, t):
        cls = type(self)
        lat, lon, alt, v_surf = body.llav(statevector(r,v), t)
        si, sj, sk = body.ijk_by_ll(lat, lon, t)
        p, atm, v_air = body.atmstate_by_lla(lat, lon, alt, t)
        return (m, r, v), (body, si, sj, sk, lat, lon, alt, v_surf), (p, atm, v-v_air), t
    
    def _do_control(self, m_rv, body_ijk_llav, pav, t):
        return self._update(m_rv, body_ijk_llav, pav, t)
    
    def _accel_dm(self, state, control):
        cls = type(self)
        ## state ##
        m_rv, body_ijk_llav, pav, t = state
        m, r, v = m_rv
        body, _, _, _, _, _, alt, _ = body_ijk_llav
        p, a, v_air = pav
        
        ## control ##
        engine_state, coefd, Tdir = control
        T, ff = cls._thrust_state(a, engine_state)
        
        accel = cls._accel_g(r, body) + cls._accel_thrust(Tdir, T, m) + cls._accel_drag(p, v_air, coefd, m)
        return accel, ff
    
    def _ode_func(self, body):
        def func(y, t):
            r = y[0:3]
            v = y[3:6]
            m = y[6]
            state = self._do_state(r,v,m,body,t)
            control = self._do_control(*state)
            a, dm = self._accel_dm(state, control)
            return concatenate((v,a,[-dm]))
        return func
    
    def _step(self, r, v, m, t, dt, odefunc):
        rvm0 = concatenate((r,v,[m]))
        rvm1 = odeint(odefunc, rvm0, list(arange(t, t+dt, 1)) + [t+dt])[-1]
        return rvm1[0:3], rvm1[3:6], rvm1[6]
    
    def _odeint(self, stv0, m0, t, dt, body):
        r0, v0 = stv0.rv
        r, v, m = self._step(r0, v0, m, t, dt, self._ode_func(body))
        return statevector(r, v), m, t + dt
    
    def _calculation_pass(self, stv, m, t, dt, body):
        return self._odeint(stv, m, t, dt, body)
    
    def _simulation_pass(self, stv, m, t, dt, body):
        self._rec = True
        res = self._odeint(stv, m, t, dt, body)
        self._rec = False
        return res


class StagedController(Controller):
    def __init__(self):
        Controller.__init__(self)
    
    def _staged_update(self, stage_name, depleted, max_twr, state):
        ## return f, Tdir, drop_stage
        ## m_rv, body_ijk_llav, pav, t = state
        raise NotImplementedError
    
    def _initialize(self, stage0):
        Controller._initialize(self)
        self._stage = stage0
        self._stage_dropped = False
    
    def _finalize(self):
        Controller._finalize(self)
        del self._stage
        del self._stage_dropped
    
    def _update(self, m_rv, body_ijk_llav, pav, t):
        stage_name, m0, m1, Tmax, isp_0, isp_1, next_stage = self._stage.statevars()
        m, _, _ = m_rv
        coefd = self._stage.coefd(m - m1)
        f, Tdir, drop_stage = self._staged_update(stage_name, (m<=m0), Tmax / (m*spconst.g), (m_rv, body_ijk_llav, pav, t))
        if stage_name == 'payload' or drop:
            self._stage_dropped = drop
            return None, coefd, zeros(3)
        else:
            self._stage_dropped = False
            return (Tmax, isp_0, isp_1, f), coefd, Tdir
    
    def sim_fullburn(self, pl, stages, stv0, body, t, dt):
        self._initialize(pl, stages)
        tf = t + dt
        while t < tf:
            name, m1, m0, Tmax, isp_0, isp_1, next = self._stage.statevars()
            if not name == 'payload':
                stage_dt = self._time_to_burnstage(stv0, body, t)
                stv0, m, _ = self._simulation_pass(stv0, m0, t, min(stage_dt, tf - t), body)
                if stage_dt > (tf - t):
                    t = tf
                    break
                else:
                    self._stage, t = next, t + stage_dt
            else: 
                stv0, m, _ = self._simulation_pass(stv0, m1, t, tf - t, body)
                break
        res = (stv0, t, self._stage.partial_depleted(m))
        self._finalize()
        return res
    
    def _time_to_burnstage(self, stv0, body, t, maxt):
        name, m1, m0, Tmax, isp_0, isp_1, _ = self._stage.statevars()
        def func(dt):
            if dt > maxt:
                return 0
            _, m, _ = self._calculation_pass(stv0, m0, t, dt, body)
            if not self.stage_dropped:
                return m - m1
            else:
                return 0
        c = smp*spconst.g/Tmax
        return brentq(func, c*min(isp_0-5, isp_1-5), c*max(isp_0+5, isp_1+5))

        
