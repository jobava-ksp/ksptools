from __future__ import print_function, division

import collections
TransferParams = collections.namedtuple('TransferParams', ['tof','cosdt','sindt','r1','r2','k','l','m','u'])
del collections

class TransferException(Exception):
    def __init__(self, msg):
        Exception.__init__(self, msg)

def getparameters(r1vec, r2vec, tof, u):
    '''Create a set of transfer parameters.'''
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
    
    return TransferParams(tof, cosdt, sindt, r1, r2, k, l, m, u)

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
    _a = min(1.0, max(-1.0, cosdE))
    dE = arccos(_a)
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
    '''Calculate orbital velocities at r1 and r2'''
    f, g, fp, gp = getfg(p, params)
    v1 = (r2-f*r1)/g
    v2 = fp*r1 + gp*v1
    return v1, v2

def toffunc(p, *params):
    '''Calculate time of flight by semi-latus rectum'''
    tof, _, _ = gaussfunc(p, TransferParams(*params))
    return tof

def toferrfunc(p, *params):
    '''Calculate difference of target time of flight and orbital time of flight.'''
    params = TransferParams(*params)
    etof = params.tof
    vtof, _ = gaussfunc(p, params)
    return etof - vtof

def toferr2func(p, *params):
    '''Calculate the square difference of target time of flight and orbital time of flight.'''
    params = TransferParams(*params)
    etof = params.tof
    vtof, _ = gaussfunc(p, params)
    return (etof - vtof)**2

def gaussfunc(p, params):
    '''Compute orbital parameters for a semi-latus rectum.'''
    from numpy import sqrt
    tof, cosdt, sindt, r1, r2, k, l, m, u = params
    f, g, fp, gp = getfg(p, params)
    a = geta(p, m, k, l)
    #e = sqrt(1-p/a)
    if a > 0:
        dE, cosdE, sindE = getdE(a, f, fp, r1, r2, u)
        tof = g + sqrt(a**3/u)*(dE-sindE)
    elif a < 0:
        dF, coshdF, sinhdF = getdF(a, f, r1)
        tof = g + sqrt((-a)**3/u)*(sinhdF-dF)
    return (tof, a)

def planar_transfer(r1vec, r2vec, tof, u, err=1e-5, maxp=1.0e+21):
    '''Calculate a transfer on a slingle plane.'''
    from numpy import array, sqrt, arcsin, empty, finfo, seterr
    from numpy.polynomial.polynomial import polyroots
    from scipy.optimize import minimize, minimize_scalar, fminbound, bisect, brentq, brenth
    
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
    elif sindt == 0.0:
        raise TransferException("True anommoly is zero or pi.")
    
    ## find points where orbit switches between hyperbolic and elliptic ##
    ## define bounds accordingly                                        ##
    bounds = [pmin] + [i for i in polyroots([-k**2, 2*k*l, 2*m-l**2]) if i >= pmin and i <= pmax] + [pmax]
    
    ## minimize over each boundry ##
    minp = []
    for i in range(len(bounds)-1):
        iter_err = err
        mn = bounds[i]
        mx = bounds[i+1]
        
        # check for small intervals #
        if mx is not None and abs(mx-mn) <= iter_err*4:
            continue
        
        # fix mn and mx for algorithm #
        if mx is None:
            mx = maxp
            p0 = min(params.r1, params.r2, maxp-iter_err)
        else:
            p0 = (mx+mn)/2
        mx -= iter_err
        mn += iter_err
        
        # check for intervals that do not intersect zero #
        mxtof, _ = gaussfunc(mx, params)
        mntof, _ = gaussfunc(mn, params)
        if (mntof - tof)*(mxtof - tof) >= 0:
            continue
        
        # bisect it down #
        p0 = brenth(toferrfunc, mn, mx, args=tuple(params), maxiter=7, disp=False)
        
        default_err = seterr(divide='raise', invalid='raise')
        optval = minimize(toferr2func, [p0], args=tuple(params), bounds=[(mn, mx)], method='TNC')
        seterr(**default_err)
        minp.append(optval)
    
    if not len(minp):
        return None
    
    ## find the very best solution ##
    p = min(minp, key=lambda opt: opt.fun).x[0]
    return getv1v2(p, r1vec, r2vec, params)

    
