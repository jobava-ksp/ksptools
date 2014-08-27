from numpy import array, matrix, identity
from numpy.linalg import norm

class RigidPosition(object):
    def __init__(self, position, orientation):
        self.position = position
        self.orientation = orientation


class RigidVelocity(object):
    def __init__(self, velocity, angular_velocity):
        self.velocity = velocity
        self.angular_velocity = angular_velocity
    
    def __mul__(self, dt):
        r = dt * self.velocity
        A = rotvec(self.axis, self.angular_velocity*dt)
        return RigidBody(r, A)
        
    def __rmul__(self, dt):
        return self.__mul__(dt)
    
    def __add__(self, other):
        v = self.velocity + other.velocity
        a = self.angular_velocity other.angular_velocity
        return RigidVelocity(v, a)
    
    def __radd__(self, other):
        return self.__add__(other)


class RigidAcceleration(object):
    def __init__(self, acceleration, angular_acceleration):
        self.acceleration = acceleration
        self.angular_acceleration = angular_acceleration
    
    def __mul__(self, dt):
        v = dt * self.acceleration
        t = dt * self.angular_acceleartion
        return RigidVelocity(v, t*self.axis)
    
    def __rmul__(self, dt):
        return self.__mul__(dt)
    
    def __add__(self, other):
        a = self.acceleration + other.acceleration
        t = self.angular_acceleartion + other.angular_acceleration
        return RigidAcceleration(a, t)
    
    def __radd__(self, other):
        return self.__add__(other)


def translate(dt, r, v, a, rp=None, vp=None):
    if rp is None:
        rp = RigidPosition(array([0,0,0]), r.orientation.copy())
    if vp is None:
        vp = RigidVelocity(array([0,0,0]), v.angular_velocit * v.axis)
    rp.position = r.position + v.velocity*dt + 0.5*(dt**2)*a.acceleration
    vp.velocity = v.velocity + a.acceleration*dt
    return rp, vp


def rotate(dt, r, v, a=None, rp=None, vp=None):
    if rp is None:
        rp = RigidPosition(r.position, matrix(identity(3)))
    if vp is None:
        vp = RigidVelocity(v.velocity, array([0,0,0]))
    
    rot = dt * v.angular_velocity + 0.5*(dt**2)*a.angular_acceleration
    rot_axis = rot/norm(rot)
    rot_t = norm(rot)
    rp.orientation = rotvec(rot_axis, rot_t) * r.orientation
    vp.angular_velocity = v.angular_velocity + dt * a.angular_acceleration
    return rp, vp


def transform(dt, r, v, a, inplace=True):
    if inplace:
        translate(dt, r, v, a, r, v)
        move(dt, r, v, a, r, v)
    else:
        r, v = translate(dt, r, v, a)
        move(dt, r, v, a, r, v)
    return r, v


class RigidBody(object):
    def __init__(self, position, velocity):
        self.position = position
        self.velocity = velocity


