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
    
    assert(tof > 0.0)
    assert(u > 0)
    
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
    import numpy as np
    assert(p > 0.0)
    while True:
        over = m*k*p
        under = ((2.*m-l**2)*p**2. + 2.*k*l*p - k**2)
        if under == 0:
            if over*under > 0:
                p += np.spacing(p)
            else:
                p -= np.spacing(p)
        else:
            break
    
    a = over/under
    return a

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
    #if dE <= 0:
    #    print(a)
    assert(dE >= 0)
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
    tof, _ = gaussfunc(p, TransferParams(*params))
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

def planar_transfer(r1vec, r2vec, tof, u, minp=0., maxp=1.0e+14):
    '''Calculate a transfer on a slingle plane.'''
    from numpy import array, sqrt, arcsin, empty, spacing, seterr
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
        bound_min = p_i
        bound_max = maxp
        pnm1 = p_i
        pn = p_i + r1
    elif sindt < 0.0:
        # dt > pi
        bound_min = minp
        bound_max = p_ii
        #p = p0 = 2*(bound_min + bound_max)/3
        #p1 = p0 + p0*1e-5
        pnm1 = minp
        pn = minp + 1e-1*minp
    elif sindt == 0.0:
        raise TransferException("True anommoly is zero or pi.")
    
    try:
        tnm1 = toffunc(pnm1, *params)
        while True:
            tn = toffunc(pn, *params)
            if abs(tn - tof) < spacing(tof):
                p = pn
                break
            pnp1 = pn + (tof - tn)*(pn - pnm1)/(tn-tnm1)
            pn, pnm1 = pnp1, pn
            tnm1 = tn
    except:
        return None
    
    '''
    ## find points where orbit switches between hyperbolic and elliptic ##
    ## define bounds accordingly                                        ##
    bounds = [bound_min] + [i for i in polyroots([-k**2, 2*k*l, 2*m-l**2]) if i > bound_min and i < bound_max] + [bound_max]
    
    ## minimize over each boundry ##
    minp = []
    for i in range(len(bounds)-1):
        mn = bounds[i]
        mx = bounds[i+1]
        
        # fix mn and mx for algorithm #
        if mx == bound_max:
            mx -= spacing((2**32)*mx)
            p0 = min((params.r1 + params.r2)/2, maxp)
        else:
            p0 = (mn+mx)/2
        mn += spacing((2**32)*mn)
        
        if mn >= mx:
            continue
        
        #default_err = seterr(divide='raise', invalid='raise')
        # bisect it down #
        #try:
        mxtof, _ = gaussfunc(mx, params)
        mntof, _ = gaussfunc(mn, params)
        #except:
        #    print("[{},{}]".format(mn,mx))
        #    print(bounds)
        #    continue
            
        if (mntof - tof)*(mxtof - tof) <= 0:
            p0 = brenth(toferrfunc, mn, mx, args=tuple(params), maxiter=9, disp=False)
        
        
        optval = minimize(toferr2func, [p0], args=tuple(params), bounds=[(mn, mx)], method='TNC')
        #seterr(**default_err)
        if optval.fun > 1:
            optval = minimize(toferr2func, [p0], args=tuple(params), bounds=[(mn, mx)], method='SLSQP')
            if optval.fun > 1:
                print("xx")
                continue
            print("{},{}".format(mn,mx))
        
        minp.append(optval)
    
    minp = [optval for optval in minp if optval.success]
    
    if not len(minp):
        return None
    
    ## find the very best solution ##
    #print(min(minp, key=lambda opt: opt.fun).fun)
    p = min(minp, key=lambda opt: opt.fun).x[0]
    '''
    return getv1v2(p, r1vec, r2vec, params)

    
