from __future__ import print_function

import collections
HohmannParams = collections.namedtuple('HohmannParams', ['tof','cosdt','sindt','r1','r2','k','l','m','u'])
del collections

def getparameters(r1vec, r2vec, tof, u):
    from ksptools.util import veccos, vecsin
    from numpy import sqrt
    from numpy.linalg import norm
    from collections import namedtuple
    
    cosdt = veccos(r1vec, r2vec)
    sindt = vecsin(r1vec, r2vec)
    r1, r2 = norm(r1vec), norm(r2vec)
    k = r1*r2*(1-cosdt)
    l = r1 + r2
    m = r1*r2*(1+cosdt)
    
    return HohmannParams(tof, cosdt, sindt, r1, r2, k, l, m, u)

def getfg(p, params):
    from numpy import sqrt, tan, arccos
    
    tof, cosdt, sindt, r1, r2, k, l, m, u = params
    tanhalfdt = tan(arccos(cosdt)/2.)
    f = 1 - r2*(1-cosdt)/p
    g = r1*r2*sindt/sqrt(u*p)
    fp = sqrt(u/p)*tanhalfdt*((1-cosdt)/p - 1/r1 - 1/r2)
    gp = 1-(r1/p)*(1-cosdt)
    return f, g, fp, gp

def geta(p, m, k, l):
    return m*k*p/((2*m-l**2)*p**2+2*k*l*p-k**2)

def getdE(a, f, fp, r1, r2):
    from numpy import sqrt, arctan2
    cosdE = 1 - r1*(1-f)/a
    sindE = r1*r2*fp/sqrt(u*a)
    dE = arctan2(cosdE,sindE)
    return dE, cosdE, sindE

def getdF(a, f, r1):
    from numpy import arccosh, sinh
    coshdF = 1 - r1*(1-f)/a
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

def toffunc(p, params):
    tof, cosdt, sindt, r1, r2, k, l, m, u = params
    f, g, fp, gp = getfg(p, params)
    
    a = geta(p, m, k, l)
    if a > 0.:
        dE, cosdE, sindE = getdE(a, f, fp, r1, r2)
        return gettofE(a, u, dE, sindE, g)
    elif a < 0.:
        dF, coshdF, sinhdF = getdF(a, f, r1)
        return gettofH(a, u, dF, sinhdF, g)

def toferrfunc(p, *params):
    params = HohmannParams(*params)
    tof = params.tof
    return (tof - toffunc(p, params))**2

def gaussorbit(r1vec, v1vec, r2vec, v2vec, tof, u):
    from numpy import array, sqrt
    from scipy.optimize import minimize, brent
    
    ## get constant parameters ##
    params = getparameters(r1vec, r2vec, tof, u)
    tof, cosdt, sindt, r1, r2, k, l, m, u = params
    
    ## pick bounds and starting semi-latus rectum ##
    if sindt > 0.0:
        # dt < pi
        pmin = k/(l+sqrt(2*m))
        pmax = float('inf')
        p0 = pmin + 1e-6
    elif sindt < 0.0:
        # dt > pi
        pmin = 1e-6
        pmax = k/(l-sqrt(2*m))
        p0 = (pmin + pmax)/2.
    
    print("(pmin, p, pmax)=({},{},{})".format(pmin,p0,pmax))
    #fmin = toferrfunc(pmin, *params)
    #fmax = toferrfunc(pmax, *params)
    #f0   = toferrfunc(p0, *params)
    #print("(fmin, f0, fmax)={}".format((fmin,f0,fmax)))
    
    
    
    #p = minimize(toferrfunc, [p0], args=tuple(params), bounds=[(pmin, pmax)], method='SLSQP')
    p = brent(toferrfunc, args=tuple(params), brack=(pmin, p0, pmax))
    f,g,fp,gp = getfg(p, params)
    vevec, vivec = getv1v2(r1vec, r2vec, f, g, fp, gp)
    return vevec-v1vec, vivec-v2vec


if __name__ == '__main__':
    import ksptools as ksp
    import ksptools.util as ksputil
    from numpy import sqrt
    from numpy.linalg import norm
    sys = ksp.loadsystem('KerbolSystem.cfg')
    s = sys['kerbol']
    k = sys['kerbin']
    e = sys['eve']
    u = s.std_g_param
    dep=13393440.0
    arv=18750672.0
    dt = arv-dep
    r1, v1 = e.getorbit().rv(dep)
    r2, v2 = k.getorbit().rv(arv)
    
    params = getparameters(r1, r2, dt, u)
    print(params)
    sindt, k, l, m = params.sindt, params.k, params.k, params.m
    print(k/(l-sqrt(2*m)))
    if sindt > 0.:
        ksputil.plotfunc(lambda n: toferrfunc(n,*params), k/(l+sqrt(2*m))+0.e-6, 3*k/(l+sqrt(2*m)))
    else:
        ksputil.plotfunc(lambda n: toferrfunc(n,*params), 0.001, k/(l-sqrt(2*m))-0.001)
    
    ve, vi = gaussorbit(r1, v1, r2, v2, arv-dep, s.std_g_param)
    print(norm(ve))
    print(norm(vi))
