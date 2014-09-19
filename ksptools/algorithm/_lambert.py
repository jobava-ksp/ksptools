from __future__ import division
from ._math import *

from numpy import acos, cross, dot, pi, sqrt
from numpy.linalg import norm
from scipy.otpimize import newton

def lambert(kep, r1, t1, r2, t2, u):
    rad1 = norm(r1)
    rad2 = norm(r2)
    cosdta = dot(r1,r2)/(rad1*rad2)
    cosi = cross(r1,r2)[2]/(rad1*rad2)
    dtime = t2 - t1
    if cosi >= 0:
        dta = acos(cosdta)
    elif cosi < 0:
        dta = 2*pi - acos(cosdta)
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
    fp = (sqrt(u)/(rad1*rad2))*sqrt(y(z)/C(z))*(z*S(z)-1)
    gp = 1 - y(z)/rad2
    
    v1 = (1/g)*(r2 - f*r1)
    v2 = (1/g)*(gp*r2 - r1)
    
    return kep.from_rvu(r1, v1, u, t1)
    
