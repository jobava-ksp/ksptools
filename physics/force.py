from numpy import cross, dot, array


class AppliedForce(object):
    def __init__(self, force, torque):
        self.torque = cross(F,r)
        self.force = (dot(F,r)*r)/dot(r,r)


def force_local(body, force, r):
    t = dot(body.A, cross(force, r))
    f = dot(body.A, (dot(force,r)*r) / dot(r,r))
    return AppliedForce(f,t)


def force_global(force, r):
    t = cross(force, r)
    f = (dot(force,r)*r)/dot(r,r)
    return AppliedForce(f,t)
    
