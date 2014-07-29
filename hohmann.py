import ksptools
import ksptools.orbit as orbit
import ksptools.util as util

def hohmann_simple(r1, r2, t1, t2):
    from numpy import cross, arccos, cos, pi
    from numpy.linalg import norm
    from scipy.integrate import romberg
    from ksptools.util import unit, reject, projmat, Ax
    
    rx = unit(r1)
    ry = unit(reject(r2,r1))
    rx = cross(rx,ry)
    plane = projmat(rx, ry, rz)
    cost = veccos(r1,r2)
    sint = vecsin(r1,r2)
    
    if sint >= 0:
        theta = arccos(cost)
    else:
        theta = 2*pi - arccos(cost)
    
    def area_func(e,a,v):
        def r(t):
            return (a*(1-e**2))/(1+e*cos(t))
        return romberg(r,v,v+dtheta)
    
    
