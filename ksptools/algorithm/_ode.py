from .._vector import statevector
from .._math import unit

from numpy import array, arange, concatenate, dot, zeros
from numpy.linalg import norm
from scipy.integrate import odeint
from scipy.optimize import brentq

import scipy.constants as spconst


def simrv(stv0, t0, t1, accel=lambda r, v, t: zeros(3)):
    rv0 = concatenate((stv0.r, stv0.v))
    def dfunc(y, t):
        r = y[0:3]
        v = y[3:6]
        a = accel(r, v, t)
        return concatenate((v, a))
    rv_array = odeint(dfunc, rv0, list(arange(t0, t1, 0.1)) + [t1])
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

    rv_array = odeint(dfunc, rv0, list(arange(t0, t1, 0.1)) + [t1])
    return statevector(rv_array[-1][0:3], rv_array[-1][3:6]), rv_array[6]


def gravity_accel_func(u):
    def gafunc(r, v, t, *args):
        return -(r/norm(r))*(u/dot(r,r))
    return gafunc


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
    
    def _odeint(self, stv, m, t, dt, body):
        r, v = stv.rv
        func = self._ode_func(body)
        r, v, m = self._step(r, v, m, t, dt, func)
        return statevector(r,v), m, t + dt
    
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
        ## return f, coefd, Tdir, drop_stage
        ## m_rv, body_ijk_llav, pav, t = state
        raise NotImplementedError
    
    def _initialize(self, payload, stages):
        Controller._initialize(self)
        self._stage = StagedController.Stage.collapse(payload, stages)
    
    def _finalize(self):
        Controller._finalize(self)
        del self._stage
    
    def _update(self, m_rv, body_ijk_llav, pav, t):
        stage_name, me, _, Tmax, isp_0, isp_1, next_stage = self._stage.statevars()
        m, _, _ = m_rv
        f, coefd, Tdir, drop = self._staged_update(stage_name, (m<=me), Tmax / (m*spconst.g), (m_rv, body_ijk_llav, pav, t))
        if stage_name == 'payload':
            return None, coefd, Tdir
        else:
            return (Tmax, isp_0, isp_1, f), coefd, Tdir
    
    def sim(self, pl, stages, stv0, body, t, dt):
        self._initialize(pl, stages)
        tf = t + dt
        while t < tf:
            name, sme, smp, Tmax, isp_0, isp_1, next = self._stage.statevars()
            if not name == 'payload':
                stage_dt = self._time_to_burnstage(stv0, body, t)
                stv0, m, _ = self._simulation_pass(stv0, sme+smp, t, min(stage_dt, tf - t), body)
                if stage_dt > (tf - t):
                    t = tf
                    break
                else:
                    self._stage, t = next, t + stage_dt
            else: 
                stv0, m, _ = self._simulation_pass(stv0, sme, t, tf - t, body)
                break
        res = (stv0, t, list(StagedController.Stage.expand_to(self._stage, m-sme)))
        self._finalize()
        return res
    
    def _time_to_burnstage(self, stv0, body, t):
        name, sme, smp, Tmax, isp_0, isp_1, _ = self._stage.statevars()
        def func(dt):
            _, m, _ = self._calculation_pass(stv0, sme + smp, t, dt, body)
            return m - sme
        c = smp*spconst.g/Tmax
        return brentq(func, c*min(isp_0-5, isp_1-5), c*max(isp_0+5, isp_1+5))
    
    @staticmethod
    def stage(name, me, mp, Tmax, isp_0, isp_1):
        return StagedController.Stage(name, me, me, mp, Tmex, isp_0, isp_1, None)
    
    @staticmethod
    def stage_by_twr(name, ms, isp_0, isp_1, Tmax, twr0, twr1):
        mp = Tmax/(spconst.g*twr0) - Tmax/(spconst.g*twr1)
        return StagedController.Stage(name, Tmax/(spconst.g*twr1), ms-mp, mp, Tmax, isp_0, isp_1, None)
    
    class Stage(object):
        def __init__(self, name, m1, me, mp, Tmax, isp_0, isp_1, next):
            self.name = name
            self.m1 = m1
            self.me = me
            self.mp = mp
            self.Tmax = Tmax
            self.isp_0 = isp_0
            self.isp_1 = isp_1
            self.next = next
        
        def statevars(self):
            return self.name, self.m1, self.mp, self.Tmax, self.isp_0, self.isp_1, self.next
        
        @staticmethod
        def collapse(pl, stages):
            macc = pl
            s0 = StagedController.Stage('payload', pl, pl, 0, 0, 0, 0, None)
            for i in range(len(stages)):
                stages[i].next = s0
                stages[i].m1 = stages[i].me + s0.m1 + s0.mp
                s0 = stages[i]
            return s0
        
        @staticmethod
        def expand(stage):
            stages = []
            while stage is not None:
                stages += [stage]
                stage = stage.next
            return stages
        
        @staticmethod
        def expand_to(stage, mp_rem):
            stage = type(stage)(stage.name, stage.m1, stage.me, mp_rem, stage.Tmax, stage.isp_0, stage.isp_1, stage.next)
            return [stage] + type(stage).expand(stage.next)

        
