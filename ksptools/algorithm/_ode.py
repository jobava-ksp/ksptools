import itertools

from .._vector import statevector
from .._math import unit

from numpy import array, arange, argmin, concatenate, diff, dot, zeros
from numpy.linalg import norm
from scipy.integrate import odeint
from scipy.optimize import brentq, minimize
from scipy.interpolate import interp1d

import scipy.constants as spconst




class State(object):
    def __init__(self, r, v, m, body, t):
        self.t = t
        self.stv = statevector(r, v)
        self.body = body
        self.llav = body.llav(self.stv, t)
        self.lat, self.lon, self.alt, self.vsurf = self.llav
        self.ijk = body.ijk_by_ll(self.lat, self.lon, t)
        self.i, self.j, self.k = self.ijk
        atm_p, atm_a, atm_v = body.atmstate_by_lla(self.lat, self.lon, self.alt, t)
        self.pav = atm_p, atm_a, v - atm_v
        self.p, self.a, self.vair = self.pav
        self.r, self.v = r, v
        self.m = m
    
    def expand(self):
        return (self.m, self.r, self.v), tuple(list(self.ijk) + list(self.llav)), self.pav, self.t


def _statevars(r, v, m, body, t):
    return State(r,v,m,body,t)


def _thruststate(a, engine_state):
    if engine_state is not None:
        Tmax, isp_0, isp_1, throttle = engine_state
        isp = isp_0 + min(1, max(0, a))*(isp_1 - isp_0)
        T = Tmax * throttle
        if isp > 0:
            return T, T/(isp * spconst.g)
        else:
            return 0, 0
    else:
        return 0, 0


class Controller(object):
    def __init__(self):
        pass
    
    def _update(self, mrv, ijkllav, pav, t):
        raise NotImplementedError
    
    def _initialize(self):
        pass
    
    def _finalize(self):
        pass
    
    def _accel_g(self, r, body):
        return -(unit(r))*(body.GM/dot(r,r))
    
    def _accel_thrust(self, Tdir, T, m):
        return (Tdir*T)/m
    
    def _accel_drag(self, p, vair, coefd, m):
        return -(0.5)*unit(vair)*p*dot(vair, vair)*coefd*8.0e-3
    
    def _accel_dm(self, body, state, control):
        engine_state, coefd, Tdir = control
        T, ff = _thruststate(state.a, engine_state)
        
        ag = self._accel_g(state.r, body)
        at = self._accel_thrust(Tdir, T, state.m)
        ad = self._accel_drag(state.p, state.vair, coefd, state.m)
        return ag + at + ad, ff
    
    def _ode_func(self, body):
        def func(y, t):
            r = y[0:3]
            v = y[3:6]
            m = y[6]
            state = State(r,v,m,body,t)
            se = state.expand()
            control = self._update(*se)
            a, dm = self._accel_dm(body, state, control)
            return concatenate((v,a,[-dm]))
        return func
    
    def _ode_initial_values(self, stv0, m0, t0, dt, step_size):
        rvm0 = concatenate((stv0.r, stv0.v, [m0]))
        x = list(arange(t0, t0+dt, step_size)) + [t0 + dt]
        return rvm0, x
    
    def _ode_split_values(self, rvm, t):
        return statevector(rvm[0:3], rvm[3:6]), rvm[6], t
    
    def _odeint(self, t0, stv0, m0, dt, body, step_size=0.1):
        rvm0, t = self._ode_initial_values(stv0, m0, t0, dt, step_size)
        rvm_array = odeint(self._ode_func(body), rvm0, t)
        return map(lambda y: self._ode_split_values(*y), zip(rvm_array,t))
    
    def _simulation_pass(self, t0, stv0, m0, dt, body, **kwargs):
        return self._odeint(t0, stv0, m0, dt, body, **kwargs)


def _minimize(data, objfunc, t0, t1, body):
    stvs, ms, ts = zip(*data)
    def yfunc(i):
        r, v = stvs[i].rv
        m = ms[i]
        t = ts[i]
        return objfunc(body, State(r,v,m,body,t))
    g0 = yfunc(1) - yfunc(0)
    g1 = yfunc(len(ts)-1) - yfunc(len(ts)-2)
    if g0 < 0 and g1 > 0:
        index = argmin([yfunc(i) for i in range(len(data))])
        r, v = stvs[index].rv
        m = ms[index]
        return True, State(r, v, m, body, ts[index])
    else:
        r, v = stvs[-1].rv
        m = ms[-1]
        return False, State(r, v, m, body, t1)


class StagedController(Controller):
    def __init__(self):
        Controller.__init__(self)
    
    def _staged_update(self, stage, depleted, max_twr, state):
        ## state = (mrv, ijkllav, pav, t)
        ## return f, Tdir, drop_stage
        raise NotImplementedError
    
    def _initialize(self, stage0):
        Controller._initialize(self)
        self._stage = stage0
        self._stage_dropped = False
        self._stage_time_dropped = None
    
    def _finalize(self):
        Controller._finalize(self)
        del self._stage
        del self._stage_dropped
        del self._stage_time_dropped
    
    def _update(self, mrv, ijkllav, pav, t):
        stage_name, m0, m1, Tmax, isp_0, isp_1, next_stage = self._stage.statevars()
        m, _, _ = mrv
        coefd = self._stage.coefd(m - m1)
        f, Tdir, drop_stage = self._staged_update(self._stage, (m<=m0), Tmax / (m*spconst.g), (mrv, ijkllav, pav, t))
        if drop_stage:
            self._stage_dropped = True
            if self._stage_time_dropped is None or t < self._stage_time_dropped:
                self._stage_time_dropped = t
            return None, coefd, zeros(3)
        elif self._stage_dropped and t >= self._stage_time_dropped:
            return None, coefd, zeros(3)
        else:
            return (Tmax, isp_0, isp_1, f), coefd, Tdir
    
    def _drop_stage(self):
        self._stage = self._stage.next
        self._stage_dropped = False
        self._stage_time_dropped = None
    
    def sim_to_objective(self, stage0, t0, stv0, body, objfunc, maxdt=900.0, step_size=0.01):
        self._initialize(stage0)
        leap_size = 1.0
        m0 = self._stage.m0
        maxt0 = maxdt + t0
        result = None
        
        while t0 < maxt0:
            data = self._simulation_pass(t0, stv0, m0, leap_size, body, step_size=step_size)
            if self._stage_dropped:
                dt = self._stage_time_dropped - t0
            else:
                dt = leap_size
            succ, state0 = _minimize(data, objfunc, t0, t0+dt, body)
            stv0, m0, t0 = state0.stv, state0.m, state0.t
            if not succ:
                if self._stage_dropped:
                    self._drop_stage()
                    m0 = self._stage.m0
            elif succ and self._stage_dropped:
                self._drop_stage()
                m0 = self._stage.m0
            else:
                result = stv0, self._stage.partial_deplete(m0), t0
                break
        self._finalize()
        return result

