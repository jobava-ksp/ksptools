
import ksptools
import ksptools.util as kspu
import ksptools.orbit as ksporbit
import hohmann as hmn
import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D

sys = ksptools.loadsystem('ToySystem.cfg')

sun = sys['sun']
akep = sys['a'].getorbit().kepler
bkep = sys['b'].getorbit().kepler
Ta = akep.period()
Tb = bkep.period()


def dohohmann(dep, dt):
    r1, v1 = akep.rv(dep)
    r2, v2 = bkep.rv(dep+dt)
    
    vh1, vh2 = hmn.gaussorbit(r1, r2, dt, sun.std_g_param)
    return ksporbit.KeplerOrbit.from_rvu(r1, vh1, sun.std_g_param, dep), r1, vh1, r2, vh2


def printahohmann(t0, dur):
    from numpy import dot, array
    from numpy.linalg import norm
    hkep, re, ve, ri, vi = dohohmann(t0,dur)

    plt.subplot(111)
    plt.plot([0],[0],'s')
    
    u = sun.std_g_param
    
    print(vars(hkep))
    print(vars(hkep.orient))
    
    r = norm(re)
    v = norm(ve)
    evec = (1./u)*((v**2-u/r)*re-dot(re,ve)*ve)
    nvec = array([1.,0.,0.])
    
    kspu.plot_semi_orbit(akep, t0, t0+dur)
    kspu.plot_semi_orbit(bkep, t0, t0+dur)
    kspu.plot_semi_orbit(hkep, t0, t0+hkep.period())
    kspu.plot_rv(re, ve)
    kspu.plot_rv(ri, vi)
    kspu.plot_rv([0,0,0],evec)
    kspu.plot_rv([0,0,0],nvec)

    plt.axis('equal')
    plt.show()
 

printahohmann(0.,3.)

