from __future__ import print_function, division

import collections
HohmannParams = collections.namedtuple('HohmannParams', ['tof','cosdt','sindt','r1','r2','k','l','m','u'])
del collections

def getparameters(r1vec, r2vec, tof, u):
    '''Create a set of gauss-hohmann parameters.'''
    from ksptools.util import veccos, vecsin, unit
    from numpy import sqrt
    from numpy.linalg import norm
    from collections import namedtuple
    
    cosdt = veccos(r1vec, r2vec)
    sindt = vecsin(r1vec, r2vec)
    
    r1, r2 = norm(r1vec), norm(r2vec)
    k = r1*r2*(1.-cosdt)
    l = r1 + r2
    m = r1*r2*(1.+cosdt)
    
    return HohmannParams(tof, cosdt, sindt, r1, r2, k, l, m, u)

def getfg(p, params):
    '''Calculate f, g, f prime, and g prime.'''
    assert(p > 0)
    assert(p < float('inf'))
    from numpy import sqrt, tan, arccos
    
    tof, cosdt, sindt, r1, r2, k, l, m, u = params
    tanhalfdt = tan(arccos(cosdt)/2.)
    f = 1. - (r2/p)*(1.-cosdt)
    g = r1*r2*sindt/sqrt(u*p)
    fp = sqrt(u/p)*tanhalfdt*((1.-cosdt)/p - 1./r1 - 1./r2)
    gp = 1. - (r1/p)*(1.-cosdt)
    return f, g, fp, gp

def geta(p, m, k, l):
    '''Calculate the semi-major axis'''
    assert(p > 0)
    return m*k*p/((2.*m-l**2)*p**2. + 2.*k*l*p - k**2)

def getdE(a, f, fp, r1, r2, u):
    '''Calculate the change in eccentric anomally for an elliptical orbit.'''
    assert(a > 0)
    from numpy import sqrt, arctan, arctan2, pi, sin, cos, arccos, arcsin
    cosdE = 1. - (r1/a)*(1.-f)
    sindE = -(r1*r2*fp)/sqrt(u*a)
    dE = arccos(cosdE)
    if sindE < 0:
        dE = 2*pi - dE
    assert(dE > 0)
    return dE, cosdE, sindE

def getdF(a, f, r1):
    '''Calculate the change in eccentric anomally for a hyperbolic orbit.'''
    assert(a < 0)
    from numpy import arccosh, sinh
    coshdF = 1. - (r1/a)*(1.-f)
    dF = arccosh(coshdF)
    sinhdF = sinh(dF)
    return dF, coshdF, sinhdF

def getv1v2(p, r1, r2, params):
    f, g, fp, gp = getfg(p, params)
    v1 = (r2-f*r1)/g
    v2 = fp*r1 + gp*v1
    return v1, v2

def toffunc(p, *params):
    tof, _, _ = gaussfunc(p, HohmannParams(*params))
    return tof

def toferrfunc(p, *params):
    params = HohmannParams(*params)
    tof = params.tof
    return (tof - toffunc(p, *params))**2

def gaussfunc(p, params):
    from numpy import sqrt
    tof, cosdt, sindt, r1, r2, k, l, m, u = params
    f, g, fp, gp = getfg(p, params)
    a = geta(p, m, k, l)
    e = sqrt(1-p/a)
    if a > 0:
        dE, cosdE, sindE = getdE(a, f, fp, r1, r2, u)
        tof = g + sqrt(a**3/u)*(dE-sindE)
    elif a < 0:
        dF, coshdF, sinhdF = getdF(a, f, r1)
        tof = g + sqrt((-a)**3/u)*(sinhdF-dF)
    return (tof, a, e)

def gaussorbit(r1vec, r2vec, tof, u):
    from numpy import array, sqrt, arcsin, empty
    from numpy.polynomial.polynomial import polyroots
    from scipy.optimize import minimize, brent, fminbound, minimize_scalar
    
    ## get constant parameters ##
    params = getparameters(r1vec, r2vec, tof, u)
    tof, cosdt, sindt, r1, r2, k, l, m, u = params
    
    ## pick bounds for the semi-latus rectum ##
    p_i  = k/(l+sqrt(2*m))
    p_ii = k/(l-sqrt(2*m))
    if sindt > 0.0:
        # dt < pi
        pmin = p_i
        pmax = None
    elif sindt < 0.0:
        # dt > pi
        pmin = 0.
        pmax = p_ii
    
    ## find points where orbit switches between hyperbolic and elliptic ##
    ## define bounds accordingly                                        ##
    bounds = [pmin] + [i for i in polyroots([-k**2, 2*k*l, 2*m-l**2]) if i >= pmin and i <= pmax] + [pmax]
    print(bounds)
    
    ## minimize over each boundry ##
    minp = []
    err = 1e-5
    for i in range(len(bounds)-1):
        mn = bounds[i]
        mx = bounds[i+1]
        if mx is None:
            p0 = 2*mn
        elif abs(mx-mn) <= err*2:
            continue
        else:
            p0 = (mn+mx)/2
            mx -= err
        optval = minimize(toferrfunc, [p0], args=tuple(params), bounds=[(mn+err,mx)], method='TNC')
        print(optval)
        minp.append(optval)
    
    ## find the very best solution ##
    p = min(minp, key=lambda opt: opt.fun).x[0]
    return getv1v2(p, r1vec, r2vec, params)

    
