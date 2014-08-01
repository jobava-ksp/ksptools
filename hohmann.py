from __future__ import print_function, division

import collections
HohmannParams = collections.namedtuple('HohmannParams', ['tof','cosdt','sindt','r1','r2','k','l','m','u'])
del collections

def getparameters(r1vec, r2vec, tof, u):
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
    from numpy import sqrt, tan, arccos
    
    tof, cosdt, sindt, r1, r2, k, l, m, u = params
    tanhalfdt = tan(arccos(cosdt)/2.)
    f = 1. - (r2/p)*(1.-cosdt)
    g = r1*r2*sindt/sqrt(u*p)
    fp = sqrt(u/p)*tanhalfdt*((1.-cosdt)/p - 1./r1 - 1./r2)
    gp = 1. - (r1/p)*(1.-cosdt)
    return f, g, fp, gp

def geta(p, m, k, l):
    return m*k*p/((2.*m-l**2)*p**2. + 2.*k*l*p - k**2)

def getdE(a, f, fp, r1, r2):
    from numpy import sqrt, arctan, arctan2, pi, sin, cos, arccos, arcsin
    cosdE = 1. - (r1/a)*(1.-f)
    sindE = -(r1*r2*fp)/sqrt(u*a)
    #dE = arctan2(cosdE,sindE)
    #dE = arctan(sindE/cosdE)
    #dE = arcsin(sindE)
    dE = arccos(cosdE)
    if sindE < 0:
        dE = 2*pi - dE
    return dE, cosdE, sindE

def getdF(a, f, r1):
    from numpy import arccosh, sinh
    coshdF = 1. - (r1/a)*(1.-f)
    dF = arccosh(coshdF)
    sinhdF = sinh(dF)
    return dF, coshdF, sinhdF

def gettofE(a, u, dE, sindE, g):
    from numpy import sqrt
    return g + sqrt(a**3/u)*(dE-sindE)

def gettofH(a, u, dF, sinhdF, g):
    from numpy import sqrt
    return g + sqrt((-a)**3/u)*(sinhdF-dF)

def getv1v2(r1, r2, f, g, fp, gp):
    v1 = (r2-f*r1)/g
    v2 = fp*r1 + gp*v1
    return v1, v2

def toffunc(p, params, returnAll=False):
    tof, cosdt, sindt, r1, r2, k, l, m, u = params
    f, g, fp, gp = getfg(p, params)
    
    a = geta(p, m, k, l)
    if a > 0.:
        dE, cosdE, sindE = getdE(a, f, fp, r1, r2)
        tof = gettofE(a, u, dE, sindE, g)
        computedE = True
    elif a < 0.:
        dF, coshdF, sinhdF = getdF(a, f, r1)
        tof = gettofH(a, u, dF, sinhdF, g)
        computedE = False
    if returnAll:
        if computedE:
            return f, g, fp, gp, a, True, dE
        else:
            return f, g, fp, gp, a, False, dF
    return tof

def toferrfunc(p, *params):
    params = HohmannParams(*params)
    tof = params.tof
    return abs(tof - toffunc(p, params))

def gaussorbit(r1vec, v1vec, r2vec, v2vec, tof, u):
    from numpy import array, sqrt, arcsin
    from scipy.optimize import minimize, brent, fminbound, minimize_scalar
    
    ## get constant parameters ##
    params = getparameters(r1vec, r2vec, tof, u)
    tof, cosdt, sindt, r1, r2, k, l, m, u = params
    
    ## pick bounds and starting semi-latus rectum ##
    p_i  = k/(l-sqrt(2*m))
    p_ii = k/(l-sqrt(2*m))
    if arcsin(sindt) > 0.0:
        # dt < pi
        pmin = p_i
        pmax = float('inf')
        p0 = pmin + 1e-6
    elif arcsin(sindt) < 0.0:
        # dt > pi
        pmin = 1e-6
        pmax = p_ii - 1e-6
        p0 = 0.7*p_ii
    
    print("(pmin, p, pmax)=({},{},{})".format(pmin,p0,pmax))
    #fmin = toferrfunc(pmin, *params)
    #fmax = toferrfunc(pmax, *params)
    #f0   = toferrfunc(p0, *params)
    #print("(fmin, f0, fmax)={}".format((fmin,f0,fmax)))
    
    
    
    #p = minimize(toferrfunc, [p0], args=tuple(params), bounds=[(pmin, pmax)], method='SLSQP')
    #p = brent(toferrfunc, args=tuple(params), brack=(pmin, pmax))
    p = fminbound(toferrfunc, pmin, pmax, args=tuple(params), disp=5)
    #p = minimize_scalar(toferrfunc, [pmin, pmax], args=tuple(params))
    f,g,fp,gp = getfg(p, params)
    vevec, vivec = getv1v2(r1vec, r2vec, f, g, fp, gp)
    return r1vec, vevec, r2vec, vivec


if __name__ == '__main__':
    import ksptools as ksp
    import ksptools.util as ksputil
    import matplotlib.pyplot as plt
    from numpy import sqrt,pi
    from numpy.linalg import norm
    sys = ksp.loadsystem('KerbolSystem.cfg')
    sun = sys['kerbol']
    kerbin = sys['kerbin']
    eve = sys['eve']
    u = sun.std_g_param
    dep=13393440.0
    arv=18750672.0
    dt = arv-dep
    ever, evev = eve.getorbit().rv(dep)
    eve_orbital_v = sqrt(eve.std_g_param/(1e+5 + eve.eq_radius))
    kerbinr, kerbinv = kerbin.getorbit().rv(arv)
    kerbin_orbital_v = sqrt(kerbin.std_g_param/(1e+5 + kerbin.eq_radius))
    
    params = getparameters(ever, kerbinr, dt, u)
    #print(params)
    sindt, k, l, m = params.sindt, params.k, params.l, params.m
    #print(k/(l-sqrt(2*m)))
    
    p_ii = k/(l-sqrt(2*m))
    left=0.1*p_ii
    right=0.9*p_ii
    
    #ksputil.plotfunc(lambda n: toffunc(n,params), k/(l+sqrt(2*m))+0.e-6, 3*k/(l+sqrt(2*m)))
    
    plt.figure()
    ax = plt.subplot(311)
    ksputil.plotfunc(lambda n: abs(toffunc(n,params) - dt), left, right, ax)
    ax = plt.subplot(312)
    for i in range(-2,3):
        ax.plot([left,right],[i*pi,i*pi])
    ksputil.plotfunc(lambda n: toffunc(n, params, True)[6], left, right, ax)
    ax = plt.subplot(313)
    ksputil.plotfunc(lambda n: toffunc(n, params, True)[4], left, right, ax)
    plt.show()
    
    _, ve, _, vi = gaussorbit(ever, evev, kerbinr, kerbinv, dt, u)
    print("from eve: {}->{}, {}dv".format(evev, ve, norm(ve-evev)))
    print("to kerbin: {}->{}, {}dv".format(kerbinv, vi, norm(vi-kerbinv)))
    
