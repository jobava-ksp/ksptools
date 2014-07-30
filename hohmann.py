from __future__ import print_function

#import ksptools
#import ksptools.orbit as orbit
#import ksptools.util as util

def hohmann_simple(r1, r2):
    from ksptools.util import arcvec, veccos
    from numpy.linalg import norm
    
    dtheta = arcvec(r1,r2)
    cost = veccos(r1,r2)
    r1_len = norm(r1)
    r2_len = norm(r2)
    rat = r1_len/r2_len
    
    if r1_len > r2_len:
        e = (rat-1)/(cost+rat)
        a = r1_len*(1-e)/(1-e**2)
    elif r1_len < r2_len:
        e = (rat-1)/(cost-rat)
        a = r1_len*(1+e)/(1-e**2)
    return e, a
    
def hohmann(r1, r2, v1, v2, dt, u):
    import numpy as np
    from ksptools.util import arcvec, veccos, Mofet, Eofet, nofretu, dMofet, dnofretu
    from numpy import pi, array
    from numpy.linalg import norm
    from scipy.optimize import minimize, fmin_tnc
    dtheta = arcvec(r1,r2)
    #print(dtheta)
    assert(dtheta > 0)
    #cosdtheta = veccos(r1,r2)
    r1_len = norm(r1)
    r2_len = norm(r2)
    
    def f(e,v):
        #if Mofet(e,v+dtheta) - Mofet(e,v) < 0.0:
        #    print("e,v,dtheta,dt: {},{},{},{}".format(e,v,dtheta,dt))
        #    print("Mv,Mv+t: {},{}".format(Mofet(e,v), Mofet(e,v+dtheta)))
        #    print("Ev,Ev+t: {},{}".format(Eofet(e,v), Eofet(e,v+dtheta)))
        #    print("nt: {}".format(nofretu(r1_len, e, v, u)*dt))
        #    assert(False)
        
        # -- full equation -- ##
        return Mofet(e,v+dtheta) - Mofet(e,v) - nofretu(r2_len, e, v+dtheta, u)*dt
        # -- mean motion equation -- ##
        '''return nofretu(r1_len, e, v, u) - nofretu(r2_len, e, v+dtheta, u)'''
    
    def fp(e,v):
        dM1e, dM1t = dMofet(e,v)
        dM2e, dM2t = dMofet(e,v+dtheta)
        _, dn1e, dn1t, _ = dnofretu(r1_len, e, v, u)
        _, dn2e, dn2t, _ = dnofretu(r2_len, e, v+dtheta, u)
        
        val = Mofet(e,v+dtheta) - Mofet(e,v) - nofretu(r2_len, e, v+dtheta, u)*dt
        
        ## -- full equation -- ##
        dval = array([dM1e - dM2e - dn2e*dt, dM1t - dM2t - dn2t*dt])
        ## -- mean motion equation -- ##
        '''dval = array([dn1e - dn2e, dn1t - dn2t])'''
        
        return dval
    
    func = lambda x: f(x[0], x[1])
    dfunc = lambda x: array(fp(x[0], x[1]))
    #cons = ({'type': 'ineq', 'func': lambda x: x[0]},
    #        {'type': 'ineq', 'func': lambda x: x[1]},
    #        {'type': 'ineq', 'func': lambda x: 2*pi - x[1]})
    
    e0 = 0.0
    t0 = 0.0
    if r1_len >= r2_len:
        e0 = r2_len/r1_len
        t0 = pi
    else:
        e0 = r1_len/r2_len
        t0 = 0.0
    
    #return f, fp, fmin_tnc(func, [e0, t0], dfunc, bounds=((0.00001, None),(0.,2*pi-0.00001)), disp=5)
    return f, fp, minimize(func, [e0, t0], method='TNC', jac=dfunc, bounds=((0.,1.0-1.0e-8),(0.,2.*pi)), options={'disp':5})
    
    #return minimize(func, [0.1, 0.], jac=dfunc, bounds=((0.,None),(0., 2*pi)))


def plot_ev(f,df,xlow=0.,xhigh=1.,ylow=None,yhigh=None):
    from numpy import array, linspace, pi, meshgrid, empty
    from scipy.optimize import minimize
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    
    if ylow is None:
        ylow = -pi
    if yhigh is None:
        yhigh = 2*pi
    
    ein = linspace(xlow, xhigh, 100)
    vin = linspace(ylow, yhigh, 100)
    
    fz = empty([len(ein), len(vin)])
    dfez = empty([len(ein), len(vin)])
    dfvz = empty([len(ein), len(vin)])
    for i in range(len(ein)):
        for j in range(len(vin)):
            fz[i,j] = f(ein[i], vin[j])
            de, dv = df(ein[i], vin[j])
            dfez[i,j] = de
            dfvz[i,j] = dv
    
    kw = {'origin':'lower', 'extent':(ylow, yhigh, xlow, xhigh)}
    
    fig = plt.figure()
    ax0 = fig.add_subplot(311)
    ax0.imshow(fz, aspect='auto', **kw)
    ax0.contour(vin,ein,fz)
    ax1 = fig.add_subplot(312)
    ax1.imshow(dfez, aspect='auto', **kw)
    ax1.contour(vin,ein,dfez)
    ax2 = fig.add_subplot(313)
    ax2.imshow(dfvz, aspect='auto', **kw)
    ax2.contour(vin,ein,dfvz)
    plt.show()

def test(r1,r2,v1,v2,dt,u):
    from numpy import pi
    f, df, mn = hohmann(r1,r2,v1,v2,dt,u)
    plot_ev(f,df,0.,1.,0.,2*pi)
    return f, df, mn
    
if __name__ == '__main__':
    from numpy import array
    import ksptools.util
    
    sys = ksptools.loadsystem('KerbolSystem.cfg')
    
    #r1 = array([1.,0.])
    #r2 = array([-1.,0.002])
    #v1 = None
    #v2 = None
    #dt = 2.5
    #u = 16.
    #_,_,m=test(r1,r2,v1,v2,dt,u)
    
    kerbol = sys['kerbol']
    kerbin = sys['kerbin']
    eve = sys['eve']
    
    r1,v1 = kerbin.orbit.rv(5793120.0)
    r2,v2 = eve.orbit.rv(8735256.0)
    dt = 8735256.0-5793120.0
    u = kerbol.std_g_param
    
    _,_,m=test(r1,r2,v1,v2,dt,u)
    print(m)
    print(r1)
    print(r2)
    print(ksptools.util.vecsin(r1,r2))
    print(ksptools.util.arcvec(r1,r2))
    
