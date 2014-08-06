import ksptools
import ksptools.util as kspu
import ksptools.orbit as ksporbit
import ksptools.transfer as ksptransfer
import ksptools.kerbaltime as ksptime

import numpy.linalg as la

#toy = ksptools.loadsystem('ToySystem.cfg')

#_sun = toy['sun']
#_akep = toy['a'].orbit.kepler
#_bkep = toy['b'].orbit.kepler
#_Ta = _akep.period()
#_Tb = _bkep.period()


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


u = kerbol.std_g_param

time_start = ksptime.KerbalDate.from_ydhms(1., 14.)
time_end = ksptime.KerbalDate.from_ydhms(1., 366)

time_start_sec = time_start.total_seconds
time_end_sec = time_end.total_seconds

ve, vi = ksptransfer.solve_transfer(u, kerbin.orbit, duna.orbit, time_start_sec, time_end_sec)
vk, vd = kerbin.orbit.v(time_start_sec), duna.orbit.v(time_end_sec)


print(vk)
print(ve)
print(vd)
print(vi)


dunar, dunav = duna.orbit.rv(time_end_sec)
kerbinr, kerbinv = kerbin.orbit.rv(time_start_sec)

solver = ksptransfer.SemiLatisSolver()
params = ksptransfer.TransferParameters(
        u,
        kerbinr,
        kerbinv,
        dunar,
        dunav,
        time_start_sec,
        time_end_sec,
        kspu.arcvec(kerbinr, dunar))

solver.start(params)

mn = solver.minp
mx = solver.maxp

#print("{}-{}".format(mn,mx))

func = lambda p: solver.tof(p, params)
#func = lambda p: (solver.tof(p, params) - (params.t1-params.t0))**2
#kspu.plotfunc(func, max(mn,1), min(mx,10e+14))
kspu.plotfunc(func, max(mn,7e+9), min(mx,1e+10))



