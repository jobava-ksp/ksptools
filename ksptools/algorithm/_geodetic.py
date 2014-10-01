from __future__ import division

from numpy import array, arccos, arcsin, cos, sin, sqrt
from numpy.linalg import norm
from scipy.optimize import minimize


def _Rlat(Re, lat, e):
    """
    :type Re: float
    :type lat: float
    :type e: float
    :rtype: float
    """
    return Re/sqrt(1-(e**2)*(sin(lat)**2))


def _gRlat(Re, lat, e):
    """
    :type Re: float
    :type lat: float
    :type e: float
    :rtype: numpy.ndarray
    """
    return array([
            1/sqrt(1-(e**2)*(sin(lat)**2)),
            Re*(e**2)*sin(lat)*cos(lat)/(1-e**2*sin(lat)**2)**(3/2),
            Re*e*sin(lat)**2/(1-e**2*sin(lat)**2)**(3/2)])


def geodetic_surface(Re, lat, e):
    """
    :type Re: float
    :type lat: float
    :type e: float
    :rtype: numpy.ndarray
    """
    rlat = _Rlat(Re, lat, e)
    return rlat*array([cos(lat), (1-(e**2))*sin(lat)])


def geodetic_latitude(r, Re, e):
    """
    :type r: numpy.ndarray
    :type Re: float
    :type e: float
    :rtype: float
    """
    r = array([norm(r[0:2]), r[2]])

    def z(lat, h):
        return (_Rlat(Re, lat, e)*(1-e**2) + h)*sin(lat)
    
    def gz(lat, h):
        gz_lat = (1-e**2)*(_Rlat(Re, lat, e)*cos(lat) + _gRlat(Re, lat, e)[1]*sin(lat)) + cos(lat)*h
        gz_h = sin(lat)
        return array([gz_lat, gz_h])
    
    def func(x):
        lat, h = x
        obj = (z(lat, h) - r[1])**2
        return obj
    
    def gfunc(x):
        lat, h = x
        g = 2*(z(lat,h) - r[1])*gz(lat, h)
        return g
        
    return minimize(func, array([arcsin(r[1]/norm(r)), norm(r) - Re]), jac=gfunc, method='BFGS').x

