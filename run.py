
import ksptools
import ksptools.util as kspu
import ksptools.orbit as ksporbit
import hohmann as hmn
import matplotlib.pyplot as plt
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


def dohohmann(akep, bkep, u, dep, dt):
    r1, v1 = akep.rv(dep)
    r2, v2 = bkep.rv(dep+dt)
    
    vh1, vh2 = hmn.gaussorbit(r1, r2, dt, u)
    return ksporbit.KeplerOrbit.from_rvu(r1, vh1, u, dep), r1, vh1, r2, vh2


def printahohmann(akep, bkep, u, t0, dur):
    from numpy import dot, array
    from numpy.linalg import norm
    hkep, re, ve, ri, vi = dohohmann(akep, bkep, u, t0,dur)
    
    r = norm(re)
    v = norm(ve)
    evec = (1./u)*((v**2-u/r)*re-dot(re,ve)*ve)
    nvec = evec
    
    fg = plt.figure()
    ax = fg.add_subplot(111, projection='3d')
    
    r1, v1 = akep.rv(t0)
    r2, v2 = bkep.rv(t0+dur)
    
    ax.plot([0],[0],[0],'s')
    kspu.plot_semi_orbit3d(akep, t0, t0+dur, ax)
    kspu.plot_semi_orbit3d(bkep, t0, t0+dur, ax)
    kspu.plot_semi_orbit3d(hkep, t0, t0+hkep.period(), ax)
    kspu.plot_rv3d(re, ve, ax)
    kspu.plot_rv3d(re, v1, ax)
    kspu.plot_rv3d(ri, vi, ax)
    kspu.plot_rv3d(ri, v2, ax)
    kspu.plot_rv3d([0,0,0],evec, ax)
    kspu.plot_rv3d([0,0,0],nvec, ax)

    ax.axis('equal')
    plt.show()
 

printahohmann(_akep, _bkep, _sun.std_g_param, 0., 11.)

