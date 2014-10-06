from __future__ import division

from .._math import C, S
from .._vector import statevector
from .._kepler import KeplerOrbit as kepler

from numpy import arccos, cos, cross, dot, pi, sin, sqrt
from numpy.linalg import norm
from scipy.optimize import newton, brentq


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
    cosz = cross(r1,r2)[2]/(rad1*rad2)
    dtime = t2 - t1
    
    if cosz >= 0:
        dta = arccos(cosdta)
    elif cosz < 0:
        dta = 2*pi - arccos(cosdta)
    
    A = sin(dta)*sqrt((rad1*rad2)/(1-cosdta))
    def y(z):
        return rad1 + rad2 + A*(z*S(z)-1)/sqrt(C(z))
    
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
            term_11 = Cz-(3/2)*(Sz/Cz)
            term_1 = (1/(2*z))*term_11 + (3/4)*((Sz**2)/Cz)
            term_2 = 3*(Sz/Cz)*sqrt(yz) + A*sqrt(Cz/yz)
            return (yz/Cz)**(3/2)*term_1 + (A/8)*term_2
    z = brentq(func, 0, 100)
    f = 1 - y(z)/rad1
    g = A * sqrt(y(z)/u)
    v1 = (1/g)*(r2 - f*r1)
    return kepler.from_statevector(statevector(r1, v1), u, t1)

