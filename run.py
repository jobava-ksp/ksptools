import datetime
import ksptools
import ksptools.util as kspu
import ksptools.orbit as ksporbit
import hohmann as hmn
import matplotlib.pyplot as plt
import pprint
from mpl_toolkits.mplot3d import Axes3D

toy = ksptools.loadsystem('ToySystem.cfg')

_sun = toy['sun']
_akep = toy['a'].getorbit().kepler
_bkep = toy['b'].getorbit().kepler
_Ta = _akep.period()
_Tb = _bkep.period()


kerbol_sys = ksptools.loadsystem('KerbolSystem.cfg')
kerbol = kerbol_sys['kerbol']
moho = kerbol_sys['moho']
eve = kerbol_sys['eve']
gilly = kerbol_sys['gilly']
kerbin = kerbol_sys['kerbin']
mun = kerbol_sys['mun']
minmus = kerbol_sys['minmus']
duna = kerbol_sys['duna']
ike = kerbol_sys['ike']
dres = kerbol_sys['dres']
jool = kerbol_sys['jool']
laythe = kerbol_sys['laythe']
vall = kerbol_sys['vall']
tyloo = kerbol_sys['tyloo']
bop = kerbol_sys['bop']
pol = kerbol_sys['pol']
eeloo = kerbol_sys['eeloo']

def dohohmann(akep, bkep, u, dep, dt):
    r1, v1 = akep.rv(dep)
    r2, v2 = bkep.rv(dep+dt)
    
    vh1vh2 = hmn.planar_transfer(r1, r2, dt, u)
    if vh1vh2:
        vh1, vh2 = vh1vh2
        return ksporbit.KeplerOrbit.from_rvu(r1, vh1, u, dep), r1, vh1, r2, vh2
    else:
        return None


def printahohmann(akep, bkep, u, t0, dur):
    from numpy import dot, array, linspace, pi
    from numpy.linalg import norm
    hkep, re, ve, ri, vi = dohohmann(akep, bkep, u, t0, dur)
    
    r = norm(re)
    v = norm(ve)
    evec = (1./u)*((v**2-u/r)*re-dot(re,ve)*ve)
    nvec = evec
    
    fg = plt.figure()
    ax = fg.add_subplot(111)
    
    r1, v1 = akep.rv(t0)
    r2, v2 = bkep.rv(t0+dur)
    
    ax.plot([0],[0],'s')
    kspu.plot_semi_orbit(akep, t0, t0+dur, ax)
    kspu.plot_semi_orbit(bkep, t0, t0+dur, ax)
    kspu.plot_semi_orbit(hkep, t0-30*dur, t0+30*dur, ax)
    #kspu.plot_rv3d(re, ve, ax)
    #kspu.plot_rv3d(re, v1, ax)
    #kspu.plot_rv3d(ri, vi, ax)
    #kspu.plot_rv3d(ri, v2, ax)
    #kspu.plot_rv3d([0,0,0],evec, ax)
    #kspu.plot_rv3d([0,0,0],nvec, ax)
    #prange = linspace(0., 1., 16., False)
    #for i in prange:
    #    r, v = hkep.rv(i*hkep.period()+t0)
    #    kspu.plot_rv(r,v,ax)
    
    pex, pey, _ = hkep.pe()
    pprint.pprint(vars(hkep))
    ax.plot([pex],[pey],'s')
    ax.axis('equal')
    plt.show()
 

def testhohmann(akep, bkep, u, t0min, t0max, dtmin, dtmax, trials):
    from itertools import product
    from numpy import linspace
    t0 = linspace(t0min, t0max, trials)
    dt = linspace(dtmin, dtmax, trials)
    cpass = 0
    cfail = 0
    for s, d in product(t0, dt):
        #try:
        r = dohohmann(akep, bkep, u, s, d)
        #except Exception as e:
        #    print("{},{}".format(s,d))
        #    print(e)
        if r is None:
            cfail += 1
        else:
            cpass += 1
    print("{}% pass".format(float(cpass)/float(cpass+cfail)))

#printahohmann(_akep, _bkep, _sun.std_g_param, 0., 0.125)

kerbin_kepler = kerbin.getorbit().kepler

n = datetime.datetime.now()
testhohmann(_akep, _bkep, _sun.std_g_param, 0., 100., 0.1, 10., 100)
print((datetime.datetime.now()-n).total_seconds())




