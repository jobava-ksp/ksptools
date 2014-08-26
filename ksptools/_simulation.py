from numpy import dot
from numpy.linalg import norm

from .util import unit
from .locallity import VectorState


class Controller(object):
    def initial_state(self, t):
        """
        :type t: float
        """
        raise NotImplementedError

    def step(self, r, v, dt):
        """
        :type r: numpy.array
        :type v: numpy.array
        :type dt: float
        """
        engines = self.getengines()
        cf = self.getcoefdrag(r, v)
        direction = self.gedirection(r, v, dt)
        return cf, direction, engines, dt

    def complete(self, r, v, mass):
        raise NotImplementedError

    @property
    def body(self):
        raise NotImplementedError


def getatmfunc(body):
    """
    :type body: ksptools.body.CelestialBody
    """
    if body.atmosphere is None:
        return lambda x: 0
    else:
        return lambda r: body.atmosphere.atm(norm(r)-body.eq_radius)


def getpfunc(body):
    """
    :type body: ksptools.body.CelestialBody
    """
    if body.atmosphere is None:
        return lambda x: 0
    else:
        return lambda r: body.atmosphere.p(norm(r)-body.eq_radius)


def accel_g(u, r):
    """
    :type u: float
    :type r: numpy.array
    """
    return -(u*unit(r))/dot(r,r)


def accel_d(cf, pfunc, r, v):
    """
    :type cf: float
    :type r: numpy.array
    :type v: numpy.array
    """
    return -unit(v)*0.5*8e-3*cf*pfunc(r)*dot(v,v)


def accel_f(d, engines, mass):
    thrust = sum(e.thrust for e in engines)
    return d * (thrust/mass)

#def physics_pass(


def run(controller, t, stepsize):
    """
    :type controller: Controller
    :type t: float
    :type stepsize: float
    """
    r, v, u, mass, body = controller.initial_state(t)
    pfunc = getpfunc(body)

    def step(dt):
        """
        :type dt: float
        """
        cf, direction, engines, dt = controller.step(r, v, dt)
        new_accel = accel_d(cf, pfunc, r, v) + accel_g(u, r) + accel_f(direction, engines, mass)
        new_v = v + new_accel*dt
        new_r = r + v*dt + 0.5*new_accel*(dt**2)
        ff = sum((e.thrust/e.isp(r)) for e in engines)
        return new_r, new_v, mass - ff*dt, t + dt

    while not controller.complete(r, v, mass):
        r, v, mass, t = step(stepsize)
    return VectorState(body, r, v, controller.body, t)