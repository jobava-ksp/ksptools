from numpy import cross, dot, array


class AppliedForce(object):
    def __init__(self, force, torque):
        self.force = force
        self.torque = torque


def force_local(body, force, r):
    t = dot(body.A, cross(force, r))
    f = dot(body.A, (dot(force,r)*r) / dot(r,r))
    return AppliedForce(f, t)


def force_global(force, r):
    t = cross(force, r)
    f = (dot(force, r)*r)/dot(r, r)
    return AppliedForce(f, t)


def torque_local(torque):
    return AppliedForce(array([0,0,0]), torque)


def torque_global(body, torque):
    t = dot(body.A, torque)
    return AppliedForce(array([0,0,0]), t)
