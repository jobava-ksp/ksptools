from itertools import count
from numpy import dot, array, cross
from numpy.linalg import norm

from .util import unit, rotvec, arcvec


class Controller(object):
    def __init__(self, vessle, phase):
        """
        :type vessle: Vessle
        """
        self.vessle = vessle
        self.refbody = vessle.state.refbody
        self.phase = phase

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
    
    def isp(self, r):
        return self.isp0 + (self.isp1 - self.isp0) * min(1, max(0, self.atm(r)))

    def accel_d(self, r, v):
        """
        :type r: numpy.ndarray
        :type v: numpy.ndarray
        """
        cf = self.dragcoef
        return -unit(v)*0.5*8e-3*cf*self.p(r)*dot(v,v)

    def accel_g(self, r):
        """
        :type r: numpy.ndarray
        """
        return -(self.refbody.std_g_param*unit(r))/dot(r,r)

    def accel_t(self, mass):
        """
        :type r: numpy.ndarray
        :type v: numpy.ndarray
        :type mass: float
        """
        return self.thrust/mass
    
    def accel_ang(self, r, v, ang, angvel):
        """
        :type r: numpy.ndarray
        :type v: numpy.ndarray
        :type ang: numpy.ndarray
        :type angvel: numpy.ndarray
        """
        destang = self.target_angle
        if destang is None:
            return array([0,0,0])
        dist = arcvec(ang, destang)
        if dist <= pi*1e-4:                 #close enough
            return array([0,0,0])
        axis = unit(cross(ang, destang))
        
        break_time = norm(angvel)/self.max_angaccel(axis)
        break_ang = norm(angvel)*break_time + 0.5*self.max_angaccel(axis)*break_time**2
        if dist <= break_ang:
            return -axis*self.max_angaccel(axis)
        elif norm(angvel) < self.max_angvel(axis):
            return axis*self.max_angaccel(axis)
        else:
            return array([0,0,0])

    def initial_state(self):
        r = self.state0.position
        v = self.state0.velocity
        t = self.state0.epoch
        ang = self.state0.orientation * uniti
        return r, v, t, ang

    def pre_step(self, r, v, ang, angvel, dt):
        """
        :type r: numpy.array
        :type v: numpy.array
        :type ang: numpy.array
        :type angvel: numpy.array
        :type dt: float
        """
        return dt, 0

    def post_step(self, r, v, ang, angvel, t):
        """
        :type r: numpy.ndarray
        :type v: numpy.ndarray
        :type ang: numpy.ndarray
        :type angvel: numpy.ndarray
        :type t: float
        """
        return True
    
    def max_angvel(self, axis):
        #TODO make it better
        return pi
    
    def max_angaccel(self, axis):
        #TODO make it better
        return pi


def run(controller, stepsize):
    """
    :type controller: Controller
    :type stepsize: float
    """
    r0, v0, t0, ang0 = controller.initial_state()
    angvel0 = array([0,0,0])

    def step(r, v, ang, angvel, t, dt):
        """
        :type r: numpy.ndarray
        :type v: numpy.ndarray
        :type ang: numpy.ndarray
        :type angvel: numpy.ndarray
        :type t: float
        :type dt: float
        """
        dt = controller.pre_step(r, v, dt)
        new_accel = controller.accel_d(r, v) + controller.accel_g(r) + ang*controller.accel_t()
        new_angaccel = controller.accel_ang(r, v, ang, angvel)
        new_v = v + new_accel*dt
        new_angvel = angvel + new_angaccel*dt
        new_r = r + v*dt + 0.5*new_accel*(dt**2)
        rvec = angvel*dt + 0.5*new_angaccel*(dt**2)
        new_ang = Ax(rotvec(unit(rvec), norm(rvec)), ang)
        return new_r, new_v, new_ang, new_angvel, t + dt

    while not controller.post_step(r0, v0, ang0, angvel0, t0):
        r0, v0, ang0, angvel0, t0 = step(r0, v0, ang0, angvel0, t0, stepsize)
    return r0, v0, t0
