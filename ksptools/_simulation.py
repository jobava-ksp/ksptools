from itertools import count
from numpy import dot, array, cross
from numpy.linalg import norm

from .util import unit, rotvec, arcvec


class Phase(object):
    def __init__(self, state, thrust=0, isp0=0, isp1=0, drymass=0, fuelmass=0, coefdrag=0, next_func=lambda x: None):
        self.state = state
        self.next = next_func
        self.thrust = thrust
        self.isp0 = isp0
        self.isp1 = isp1
        self.drymass = drymass
        self.fuelmass = fuelmass
        self.coefdrag = coefdrag
        self.target_angle = None
    
    def update(self, r, v, ang, angvel, t, mass):
        return self


class BurnPhase(Phase):
    def __init__(self, thrust, isp0, isp1, drymass, fuelmass, coefdrag, next_func=lambda x: None):
        Phase.__init__(self, 'burn', thrust, isp0, isp1, drymass, fuelmass, coefdrag, next_func)
    
    def drawfuel(self, dmass):
        self.fuelmass = max(0, self.fuelmass - dmass)
        if dmass > self.fuelmass:
            return self.fuelmass
        else:
            return dmass
    
    def update(self, r, v, ang, angvel, t, mass):
        if self.fuelmass <= 0:
            return self.next(self)
        else:
            return self


class PreLaunchPhase(Phase):
    def __init__(self, state, burn_phase):
        Phase.__init__(self, 'prelaunch', next_func = lambda x: burn_phase)
    
    def direction(self, r, v, t):
        return unit(r)


class Controller(object):
    def __init__(self, state, phase):
        """
        :type state: ksptools.locallity.State
        """
        self.state0 = state
        self.refbody = state.refbody
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
        mass = self.phase.drymass + self.phase.fuelmass
        direction = self.phase.direction(r,v,t)
        return r, v, t, direction, mass

    def pre_step(self, r, v, ang, angvel, dt):
        """
        :type r: numpy.array
        :type v: numpy.array
        :type ang: numpy.array
        :type angvel: numpy.array
        :type dt: float
        """
        if self.phase.state == 'burn':
            ff = self.isp(r) /self.thrust
            max_dmass = dt*ff
            avb_dmass = self.phase.drawfuel(max_dmass)
            if max_dmass >= avb_dmass:
                self.phase = self.phase.next()
                dt, dmass = avb_dmass/ff, -avb_dmass
            else:
                dt, dmass = dt, -max_dmass
        return dt, dmass

    def post_step(self, r, v, ang, angvel, t, mass):
        self.phase = self.phase.update(r, v, ang, angvel, t, mass)
        if self.phase is None:
            return True
        self.mass = self.phase.fuelmass + self.phase.drymass
        self.thrust = self.phase.thrust
        self.isp0 = self.phase.isp0
        self.isp1 = self.phase.isp1
        self.coefdrag = self.phase.coefdrag
        self.target_angle = self.phase.target_angle
        return False
    
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
    r0, v0, t0, mass0, ang0 = controller.initial_state()
    angvel0 = array([0,0,0])

    def step(r, v, ang, angvel, t, dt, mass):
        """
        :type r: numpy.ndarray
        :type v: numpy.ndarray
        :type ang: numpy.ndarray
        :type angvel: numpy.ndarray
        :type t: float
        :type dt: float
        :type mass: float
        """
        dt, dmass = controller.pre_step(r, v, dt)
        new_accel = controller.accel_d(r, v) + controller.accel_g(r) + ang*controller.accel_t(mass)
        new_angaccel = controller.accel_ang(r, v, ang, angvel)
        new_v = v + new_accel*dt
        new_angvel = angvel + new_angaccel*dt
        new_r = r + v*dt + 0.5*new_accel*(dt**2)
        rvec = angvel*dt + 0.5*new_angaccel*(dt**2)
        new_ang = Ax(rotvec(unit(rvec), norm(rvec)), ang)
        return new_r, new_v, new_ang, new_angvel, mass + dmass, t + dt

    while not controller.post_step(r0, v0, ang0, angvel0, t0, mass0):
        r0, v0, ang0, angvel0, mass0, t0 = step(r0, v0, ang0, angvel0, t0, stepsize, mass0)
    return r0, v0, t0
