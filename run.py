
import ksptools
import ksptools.util as kspu
import matplotlib.pyplot as plt

sys = ksptools.loadsystem('ToySystem.cfg')

a = sys['a'].getorbit().kepler
b = sys['b'].getorbit().kepler
Ta = a.period()
Tb = b.period()

kspu.plot_semi_orbit(a,0,Ta)
kspu.plot_semi_orbit(b,0,Tb)
kspu.plot_rv(a,0.125*Ta)
kspu.plot_rv(b,0.125*Tb)
kspu.plot_rv(a,0.375*Ta)
kspu.plot_rv(b,0.375*Tb)
kspu.plot_rv(a,0.625*Ta)
kspu.plot_rv(b,0.625*Tb)
kspu.plot_rv(a,0.875*Ta)
kspu.plot_rv(b,0.875*Tb)
plt.plot([0],[0],'s')
plt.axis('equal')
plt.show()

