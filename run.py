from __future__ import division

import ksptools as ksp
import ksptools.orbit as kspo
import ksptools.part as ksp_part
import ksptools.util as kspu
import ksptools.kerbaltime as kspt
import ksptools.transfer as ksptransfer

import numpy as np
import numpy.linalg as npla

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from kerbolsystem import *

#return vp, rpe, rsoi, eta - delta_time
def plot(sol, vsoi, tpe, u):
    vpe, rpe, rsoi, tsoi = sol
    kep = kspo.KeplerOrbit.from_rvu(rsoi, vsoi, u, tsoi)
    
    start_time = min(tpe, tsoi)
    end_time = max(tpe, tsoi)
    dur = max(tpe, tsoi) - start_time
    r0 = kep.r(start_time)
    r1 = kep.r(start_time - dur)
    r2 = kep.r(start_time + dur)
    
    xp, yp, _ = rpe
    x0, y0, _ = r0
    x1, y1, _ = r1
    x2, y2, _ = r2
    
    plt.plot([0],[0],'*')
    plt.plot([xp],[yp],'sg')
    plt.plot([x0],[y0],'sk')
    #plt.plot([x1],[y1],'sr')
    plt.plot([x2],[y2],'sr')
    
    kspu.plot_semi_orbit(kep, start_time, start_time + dur, plt)
    
    print("vp: {}".format(vpe))
    print("rp: {}, {}, ({})".format(npla.norm(rpe), npla.norm(rpe)-kerbin.eq_radius, rpe))
    print("rsoi: {}, ({})".format(npla.norm(rsoi), rsoi))
    print("tsoi <-> tpe: {} <-> {}".format(tsoi, tpe))
    

def insert_eject(vsoi):
    t1 = kspt.toseconds(4)
    eject_sol = ksptransfer.solve_flyby(kerbin.std_g_param, 8.0e+4 + kerbin.eq_radius, vsoi, kerbin.soi, t1, insertion=False)
    insert_sol = ksptransfer.solve_flyby(kerbin.std_g_param, 8.0e+4 + kerbin.eq_radius, vsoi, kerbin.soi, t1, insertion=True)
    
    plot(eject_sol, vsoi, t1, kerbin.std_g_param)
    plot(insert_sol, vsoi, t1, kerbin.std_g_param)
    plt.axis('equal')
    plt.show()

#insert_eject(np.array([400,0,-40]))

def make_parts():
    liquid_fuel = ksp_part.ResourceType('liquid_fuel', 'Liquid Fuel', 5, 0.8, 'stage')
    oxidizer = ksp_part.ResourceType('oxidizer', 'Oxidizer', 5, 0.18, 'stage')
    monoprop = ksp_part.ResourceType('monoprop', 'Mono Propellant', 4, 1.2, 'universal')

    #CrewPodPartType(name, title, radialsize, cost, mass, ceofdrag, kerbals, torque, resources)
    mk1 = ksp_part.CrewPodPartType(
        'mk1pod', 'Command Pod Mk-1', 'small', 588, 0.8e+3, 0.2, 1, 5.0e+3, [(10, monoprop)])

    #CutePartType(name, title, redialsize, cost, mass, coefdrag, semidrag, fulldrag)
    mk16chute = ksp_part.ChutePartType(
        'mk16chute', 'Mk16 Parachute', 'tiny', 422, 0.1e+3, 0.22, 1.0, 500.0, 0.01, 500)

    #DecouplerPartType(name, title, radialsize, cost, mass, coefdrag, ejection_momentum)
    tr18a = ksp_part.DecouplerPartType(
        'tr18a', 'TR-18A Stack Decoupler', 'small', 400, 0.5e+3, 0.2, 250.0e+3)

    #FuelTankPartType(name, title, radialsize, cost, mass, coefdrag, resources)
    flt100 = ksp_part.FuelTankPartType(
        'flt100', 'FL-T100 Fuel Tank', 'small', 204.1, 0.06e+3, 0.2, [(45, liquid_fuel), (55, oxidizer)])

    #EnginePartType(name, title, radialsize, cost, mass, coefdrag, fueltypes, minthrust, maxthrust, isp, ispatm)
    lvt45 = ksp_part.EnginePartType(
        'lvt45', 'LV-T45 Liquid Fule Engine', 'small', 950, 1.5e+3, 0.2, [(0.45, liquid_fuel), (0.55, oxidizer)], 0, 200.0e+3, 370, 320)
    
make_parts()
