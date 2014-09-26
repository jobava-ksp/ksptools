from __future__ import division

from .._math import *
from .._vector import statevector
from .._kepler import KeplerOrbit as kepler

from numpy import arccos, cross, dot, pi, sqrt
from numpy.linalg import norm
from scipy.optimize import newton


def lambert(r1, t1, r2, t2, u):
    """
    :type r1: numpy.ndarray
    :type t1: float
    :type r2: numpy.ndarray
    :type t2: float
    :type u: float
    :rtype: ksptools._kepler.KeplerOrbit
    """
    rad1 = norm(r1)
    rad2 = norm(r2)
    cosdta = dot(r1,r2)/(rad1*rad2)
    cosi = cross(r1,r2)[2]/(rad1*rad2)
    dtime = t2 - t1
    if cosi >= 0:
        dta = arccos(cosdta)
    elif cosi < 0:
        dta = 2*pi - arccos(cosdta)
    else:
        raise NotImplementedError
    
    A = sin(dta)*sqrt((rad1*rad2)/(1-cosdta))
    y = lambda z: rad1 + rad2 + A*(z*S(z)-1)/sqrt(C(z))
    
    def func(z):
        yz = y(z)
        return (yz/C(z))**(3/2)*S(z) + A*sqrt(yz) - sqrt(u)*dtime
    
    def funcp(z):
        if z == 0:
            y0 = y(0)
            return (sqrt(2)/40)*y0**(3/2)+(A/8)*(sqrt(y0) + A*sqrt(1/(2*y0)))
        else:
            yz = y(z)
            Cz = C(z)
            Sz = S(z)
            lterm = (yz/Cz)**(3/2)*((1/(2*z))*(Cz-(3/2)*(Sz/Cz)) + (3/4)*(Sz**2/Cz))
            rterm = (A/8)*(3*(Sz/Cz)*sqrt(yz) + A*sqrt(Cz/yz))
            return lterm + rterm
    
    z = newton(func, 0, funcp)
    f = 1 - y(z)/rad1
    g = A * sqrt(y(z)/u)
    #fp = (sqrt(u)/(rad1*rad2))*sqrt(y(z)/C(z))*(z*S(z)-1)
    #gp = 1 - y(z)/rad2
    
    v1 = (1/g)*(r2 - f*r1)
    #v2 = (1/g)*(gp*r2 - r1)
    
    return kepler.from_statevector(statevector(r1, v1), u, t1)

