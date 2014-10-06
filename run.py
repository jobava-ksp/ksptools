from ksptools import *



ksys = kerbolsystem
t0 = datetosec('Y1 D1, 00:00:00')
t1 = datetosec('Y2 D72, 3:48:00')

print(t0)
print(t1)

from ksptools.algorithm._manuver import solve_ejection_prograde, solve_insertion_prograde

def direct_transfer(bodyA, rpe0, t0, bodyB, rpe3, t3, sun):
    stvi = bodyA.statevector(t0)
    stvf = bodyB.statevector(t3)
    kepT = kepler.lambert(stvi.r, t0, stvf.r, t3, sun.GM)
    stv0 = kepT.statevector_by_time(t0)
    stv3 = kepT.statevector_by_time(t3)
    kepI, kepA, t0, t1 = solve_ejection_prograde(bodyA, stv0, t0, rpe0, sun)
    kepF, kepB, t2, t3 = solve_insertion_prograde(bodyB, stv3, t3, rpe3, sun)
    stv1 = kepT.statevector_by_time(t1)
    stv2 = kepT.statevector_by_time(t2)
    kepT = kepler.lambert(stv1.r, t1, stv2.r, t2, sun.GM)
    dv_eject = kepA.statevector_by_time(t0).v - kepI.statevector_by_time(t0).v
    dv_insert = kepF.statevector_by_time(t3).v - kepB.statevector_by_time(t3).v
    return (dv_eject, dv_insert), [(kepA, t0, t1), (kepT, t1, t2), (kepB, t2, t3)]

import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt
from scipy.optimize import minimize

def porkchop(bodyA, rpe0, bodyB, rpe3, ts, sun, res):
    dept_range = min(bodyA.synodic(bodyB), bodyA.period)
    rap = max(bodyA.rap, bodyB.rap)
    rpe = min(bodyA.rpe, bodyA.rpe)
    Th = 2*(np.pi/np.sqrt(sun.GM))*((rap+rpe)/2)**(3.0/2.0)
    dt_min = max(Th/2 - bodyB.period, Th/4)
    dt_max = dt_min + min(2*bodyB.period, Th/2)
    
    dept = np.linspace(ts, ts + dept_range, res)
    durt = np.linspace(ts + dt_min, ts + dt_max, res)
    
    xv, yv = np.meshgrid(dept, durt)
    z = np.empty((res,res))
    
    def dvfunc(x):
        t0, dur = x
        #try:
        #    (dve, dvi), _ = direct_transfer(bodyA, rpe0, t0, bodyB, rpe3, t0+dur, sun)
        #except Exception as e:
        #    print('{} - {}'.format(sectodate(t0), sectotime(dur)))
        #    print(str(e))
        #    exit()
            
        (dve, dvi), _ = direct_transfer(bodyA, rpe0, t0, bodyB, rpe3, t0+dur, sun)
        return la.norm(dve) + la.norm(dvi)
    
    for i in range(len(dept)):
        for j in range(len(durt)):
            z[i,j] = dvfunc((dept[i], durt[j]))
    #print(minimize(dvfunc, [ts + dept_range/2, ts + (dt_min + dt_max)/2]))
    #print(z[0])
    print('\n'.join(''.join(map(lambda x: 'X' if x == 0 else ' ', z[i])) + '|' for i in range(res)))

porkchop(ksys['kerbin'], 6.85e+5, ksys['duna'], 3.7e+5, t0, ksys['kerbol'])
#(dve, dvi), _ = direct_transfer(ksys['kerbin'], 6.85e+5, t0, ksys['duna'], 3.7e+5, t1, ksys['kerbol'])
#print('{} -> {}'.format(la.norm(dve), la.norm(dvi)))

def x(arr):
    return [arr[i].r[0] for i in range(len(arr))]

def y(arr):
    return [arr[i].r[1] for i in range(len(arr))]

def z(arr):
    return [arr[i].r[2] for i in range(len(arr))]



