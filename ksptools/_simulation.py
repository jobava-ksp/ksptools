from numpy import dot
from numpy.linalg import norm

from .util import unit


class SimEngine(object):
    def __init__(self, thrust, isp, isp_atm, requestfuel_func, drawfuel_func, deplete_func):
        self.maxthrust = thrust
        self.throttle = 0
        self.isp_0 = isp
        self.isp_1 = isp_atm
        self.requestfuel = requestfuel_func
        self.drawfuel = drawfuel_func
        self.deplete = deplete_func

    @property
    def thrust(self):
        return self.thrust * self.throttle

    def isp(self, atm):
        return self.isp_0 - min(1, max(0, atm))*(self.isp_0 - self.isp_1)

    def ff(self, atm):
        return self.isp(atm)/self.thrust

    def request_time(self, dt, atm):
        req_fuel = dt*self.ff(atm)
        max_fuel = self.requestfuel(req_fuel)
        return max(0, dt*(max_fuel/req_fuel))

    def draw_time(self, dt, atm):
        if dt <= 0:
            self.deplete()
        else:
            draw_fuel = dt*self.ff(atm)
            self.drawfuel(draw_fuel)


class Controller(object):
    def __init__(self, state):
        """
        :type state: ksptools.locallity.State
        """
        self.state0 = state
        self.refbody = state.refbody

    def p(self, r):
        if self.refbody.atmosphere is None:
            return 0
        else:
            return self.refbody.atmosphere.p(norm(r)-self.refbody.eq_radius)

    def atm(self, r):
        if self.refbody.atmosphere is None:
            return 0
        else:
            return self.refbody.atmosphere.atm(norm(r)-self.refbody.eq_radius)

    def accel_d(self, r, v):
        """
        :type r: numpy.array
        :type v: numpy.array
        """
        cf = self.get_dragcoef(r, v)
        return -unit(v)*0.5*8e-3*cf*self.p(r)*dot(v,v)

    def accel_g(self, r):
        """
        :type r: numpy.array
        """
        return -(self.refbody.std_g_param*unit(r))/dot(r,r)

    def accel_t(self, r, v, dt, mass):
        """
        :type r: numpy.array
        :type v: numpy.array
        :type mass: float
        """
        thrust = sum(e.thrust for e in self.get_engines())
        return self.get_direction(r, v, dt) * (thrust/mass)

    def initial_state(self):
        r = self.state0.position
        v = self.state0.velocity
        t = self.state0.epoch
        mass = self.get_initialmass()
        return r, v, t, mass

    def step(self, r, v, dt):
        """
        :type r: numpy.array
        :type v: numpy.array
        :type dt: float
        """
        atm = self.atm(r)
        dt = min(e.request_time(dt, atm) for e in self.get_engines())
        dmass = -sum(e.ff(atm)*dt for e in self.get_engines())
        for e in self.get_engines():
            e.draw_time(dt)
        return dt, dmass

    def logic_step(self, r, v, t):
        raise NotImplementedError

    def get_engines(self):
        raise NotImplementedError

    def get_initialmass(self):
        raise NotImplementedError

    def get_direction(self, r, v, dt):
        raise NotImplementedError

    def get_dragcoef(self, r, v):
        raise NotImplementedError


def run(controller, stepsize):
    """
    :type controller: Controller
    :type stepsize: float
    """
    r0, v0, t0, mass0 = controller.initial_state()

    def step(r, v, t, mass, dt):
        """
        :type dt: float
        """
        dt, dmass = controller.step(r, v, dt)
        new_accel = controller.accel_d(r, v) + controller.accel_g(r) + controller.accel_t(r, v, dt, mass)
        new_v = v + new_accel*dt
        new_r = r + v*dt + 0.5*new_accel*(dt**2)
        return new_r, new_v, mass + dmass, t + dt

    while not controller.logic_step(r0, v0, t0):
        r0, v0, mass0, t0 = step(r0, v0, t0, mass0, stepsize)
    return r0, v0, t0