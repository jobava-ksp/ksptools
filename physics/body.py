from .node import TreeNode
from numpy import array, mat, asmatrix, identity, dot, cos, sin
from numpy.linalg import norm


def _rotmat(u):
    """
    :type u: numpy.ndarray
    """
    t = norm(u)
    x, y, z = (u/t)
    ct, st = cos(t), sin(t)
    return mat([[ct+x*x*(1-ct),   x*y*(1-ct)-z*st, x*z*(1-ct)+y*st],
                [y*x*(1-ct)+z*st, ct+y*y*(1-ct),   y*z*(1-ct)-x*st],
                [z*x*(1-ct)-y*st, z*y*(1-ct)+x*st, ct+z*z*(1-ct)  ]])


class RigidPosition(object):
    def __init__(self, position, orientation):
        """
        :type position: numpy.ndarray
        :type orientation: numpy.ndarray
        """
        self.position = position
        self.orientation = orientation
    
    @staticmethod
    def translate(dt, r, v, a, rp=None, vp=None):
        """
        :type dt: float
        :type r: RigidPosition
        :type v: RigidVelocity
        :type a: RigidAcceleration
        :type rp: RigidPosition
        :type vp: RigidVelocity
        """
        if rp is None:
            rp = RigidPosition(array([0,0,0]), r.orientation.copy())
        if vp is None:
            vp = RigidVelocity(array([0,0,0]), v.angular_velocit * v.axis)
        rp.position = r.position + v.velocity*dt + 0.5*(dt**2)*a.acceleration
        vp.velocity = v.velocity + a.acceleration*dt
        return rp, vp
    
    @staticmethod
    def rotate(dt, r, v, a=None, rp=None, vp=None):
        """
        :type dt: float
        :type r: RigidPosition
        :type v: RigidVelocity
        :type a: RigidAcceleration
        :type rp: RigidPosition
        :type vp: RigidVelocity
        """
        if rp is None:
            rp = RigidPosition(r.position, asmatrix(identity(3)))
        if vp is None:
            vp = RigidVelocity(v.velocity, array([0,0,0]))
        
        rot = dt * v.angular_velocity + 0.5*(dt**2)*a.angular_acceleration
        rp.orientation = _rotmat(rot) * r.orientation
        vp.angular_velocity = v.angular_velocity + dt * a.angular_acceleration
        return rp, vp

    @staticmethod
    def transform(dt, r, v, a, inplace=True):
        """
        :type dt: float
        :type r: RigidPosition
        :type v: RigidVelocity
        :type a: RigidAcceleration
        :type inplace: bool
        """
        if inplace:
            translate(dt, r, v, a, r, v)
            rotate(dt, r, v, a, r, v)
        else:
            r, v = translate(dt, r, v, a)
            rotate(dt, r, v, a, r, v)
        return r, v


translate = RigidPosition.translate
rotate = RigidPosition.rotate
transform = RigidPosition.transform


class RigidVelocity(object):
    def __init__(self, velocity, angular_velocity):
        """
        :type velocity: numpy.ndarray
        :type angular_velocity: numpy.ndarray
        """
        self.velocity = velocity
        self.angular_velocity = angular_velocity
    
    def __mul__(self, dt):
        """
        :type dt: float
        """
        r = dt * self.velocity
        A = _rotmat(r)
        return RigidPosition(r, A)
        
    def __rmul__(self, dt):
        """
        :type dt: float
        """
        return self.__mul__(dt)
    
    def __add__(self, other):
        """
        :type other: RigidVelocity
        """
        v = self.velocity + other.velocity
        a = self.angular_velocity + other.angular_velocity
        return RigidVelocity(v, a)
    
    def __radd__(self, other):
        """
        :type other: RigidVelocity
        """
        return self.__add__(other)


class RigidAcceleration(object):
    def __init__(self, acceleration, angular_acceleration):
        """
        :type acceleration: numpy.ndarray
        :type angular_acceleration: numpy.ndarray
        """
        self.acceleration = acceleration
        self.angular_acceleration = angular_acceleration
    
    def __mul__(self, dt):
        """
        :type dt: float
        """
        v = dt * self.acceleration
        w = dt * self.angular_acceleration
        return RigidVelocity(v, w)
    
    def __rmul__(self, dt):
        """
        :type dt: float
        """
        return self.__mul__(dt)
    
    def __add__(self, other):
        """
        :type other: RigidAcceleration
        """
        a = self.acceleration + other.acceleration
        t = self.angular_acceleration + other.angular_acceleration
        return RigidAcceleration(a, t)
    
    def __radd__(self, other):
        """
        :type other: RigidAcceleration
        """
        return self.__add__(other)


class RigidBody(TreeNode):
    def __init__(self, position, velocity, mass):
        """
        :type position: RigidPosition
        :type velocity: RigidVelocity
        :type mass: float
        """
        TreeNode.__init__(self)
        self.mass = mass
        self._position = position
        self._velocity = velocity
    
    def _get_x(self):
        return self._position.position
    
    def _set_x(self, x):
        self._position.position = x
    
    def _get_v(self):
        return self._velocity.velocity
    
    def _set_v(self, v):
        self._velocity.velocity = v
    
    def _get_A(self):
        return self._position.orientation
    
    def _set_A(self, A):
        self._position.orientation = A
    
    def _get_w(self):
        return self._velocity.angular_velocity
    
    def _set_w(self, w):
        self._velocity.angular_velocity = w
    
    def _coefdrag_local(self, d):
        """
        :type d: numpy.ndarray
        """
        raise NotImplementedError
    
    def _area_local(self, d):
        """
        :type d: numpy.ndarray
        """
        raise NotImplementedError
    
    def _inertia_local(self, axis):
        """
        :type axis: numpy.ndarray
        """
        raise NotImplementedError
    
    def transform_newton(self, dt, accel, angular_accel):
        """
        :type dt: float
        :type accel: numpy.ndarray
        :type angular_accel: numpy.ndarray
        """
        a = RigidAcceleration(accel, angular_accel)
        translate(dt, self._position, self._velocity, a)
    
    def coefdrag(self, d):
        """
        :type d: numpy.ndarray
        """
        d = dot(self.A.T, d)
        return self._coefdrag_local(d)
    
    def area(self, d):
        """
        :type d: numpy.ndarray
        """
        d = dot(self.A.T, d)
        return self._area_local(d)
    
    def inertia(self, axis):
        """
        :type axis: numpy.ndarray
        """
        axis = dot(self.A.T, axis)
        return self._inertia_local(axis)
    
    x = property(_get_x, _set_x)
    v = property(_get_v, _set_v)
    A = property(_get_A, _set_A)
    w = property(_get_w, _set_w)
    

