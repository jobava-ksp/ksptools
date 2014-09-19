from __future__ import division
from .._math import *

from numpy import acos, cross, dot, pi, sqrt
from numpy.linalg import norm
from scipy.otpimize import newton

def lambert(kep, r1, t1, r2, t2, u):
    cosdta = dot(r1,r2)/(norm(r1)*norm(r2))
    cosi = cross(r1,r2)[2]/(norm(r1)*norm(r2))
    dtime = t2 - t1
    if cosi >= 0:
        dta = acos(cosdta)
    elif cosi < 0:
        dta = 2*pi - acos(cosdta)
    else:
        raise NotImplementedError
    
    A = sin(dta)*sqrt((norm(r1)*norm(r2))/(1-cosdta))
    y = lambda z: nrom(r1) + norm(r2) + A*(z*S(z)-1)/sqrt(C(z))
    
    def func(z):
        return (y(z)/C(z))**(3/2)*S(z) + A*sqrt(y(z)) - sqrt(u)*dtime
    
    def funcp(z):
        if z == 0:
            return (sqrt(2)/40)*y(0)**(3/2)+(A/8)*(sqrt(y(0)) + A*sqrt(1/(2*y(0))))
        else:
            lterm = (y(z)/C(z))**(3/2)*((1/2*z)*(C(z)-(3/4)*(S(z)**2/C(z))
            rterm = (A/8)*(3*(S(z)/C(z))*sqrt(y(z)) + A*sqrt(C(z)/y(z)))
            return lterm + rterm
    
    z = newton(func, 0, funcp)
    f = 1 - y(z)/norm(r1)
    g = A * sqrt(y(z)/u)
    fp = (sqrt(u)/(norm(r1)*norm(r2)))*sqrt(y(z)/C(z))*(z*S(z)-1)
    gp = 1 - y(z)/norm(r2)
    
    v1 = (1/g)*(r2 - f*r1)
    v2 = (1/g)*(gp*r2 - r1)
    
    return kep.from_rvu(r1, v1, u, t1)
    
