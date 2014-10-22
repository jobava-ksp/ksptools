from __future__ import division
import itertools

from . import _ode
from .._kepler import KeplerOrbit as kepler
from .._vector import statevector
from .._math import rotaxis, unit, unitk
#from ..path._path import Path as path

from numpy import arccos, array, cross, dot, imag, pi, roots, sin, sqrt
from numpy.linalg import norm
from scipy.optimize import minimize


def phasechange(r0, r1, time_anomaly, n, mu):
    dt = time_anomaly/n
    period = ((2*pi)/sqrt(mu))*((r0+r1)/2)**(3/2)
    rp = 2*((period+dt)*sqrt(mu)/(2*pi))**(2/3) - r0
    hs = sqrt(mu*(r0*r1)/(r0+r1))
    hp = sqrt(mu*(r0*rp)/(r0+rp))
    return rp, (hp - hs)/r0, n*(period+dt)


def incline(r0, delti, mu):
    def dv(rp):
        ht = sqrt(mu*(r0*rp)/(r0+rp))
        hc = sqrt(0.5*mu*r0)
        return 2*(abs(ht-hc)/r0) + 2*(ht/rp)*sin(delti/2)
    def dvp(rp):
        return (dv(rp + 50) - dv(rp-50))/100
    res = minimize(dv, r0, jac=dvp)
    #print(res)
    rp = res.x[0]
    ht = sqrt(mu*(r0*rp)/(r0+rp))
    hc = sqrt(0.5*mu*r0)
    return rp, (ht-hc)/r0, 2*(ht/rp)*sin(delti/2), dv(rp)


def _flyby_parameters(body, stv0, t0, rpe, sun):
    ### orbital parameters (perifocal) ###
    #vesc = sun.relative_statevector(body, stv0, t0).v ## TODO: fix this
    vesc = stv0.v - body.statevector(t0).v
    e = 1 + rpe * dot(vesc, vesc)/body.GM
    vpe = sqrt(dot(vesc,vesc) + 2*body.GM/rpe)
    vce = sqrt(body.GM/rpe)
    
    ### orientation ###
    B = arccos(1/e)            # angle of assymptote
    normal = unit(cross(cross(vesc, unitk), vesc))

    return vesc, vpe, vce, B, normal


def solve_ejection_prograde(body, stv0, t0, rpe, sun=None):
    if sun is None:
        sun = body.parent_node
    vesc, vpe, vce, B, normal = _flyby_parameters(body, stv0, t0, rpe, sun)

    r0 = rpe * dot(rotaxis(normal, B), (unit(-vesc))).A1
    v0 = vpe * unit(cross(normal, r0))
    vi = vce * unit(cross(normal, r0))
    
    kepO = kepler.from_statevector(statevector(r0, vi), body.GM, t0)
    kepE = kepler.from_statevector(statevector(r0, v0), body.GM, t0)
    taf = kepE.true_anomaly_by_distance(body.soi)
    return kepO, kepE, t0, t0 + kepE.time_anomaly_by_ta(taf)


def solve_insertion_prograde(body, stv0, t0, rpe, sun=None):
    if sun is None:
        sun = body.parent_node
    vesc, vpe, vce, B, normal = _flyby_parameters(body, stv0, t0, rpe, sun)

    r0 = rpe * dot(rotaxis(normal, -B), (-unit(vesc))).A1
    v0 = vpe * unit(cross(normal, r0))
    vi = vce * unit(cross(normal, r0))
    
    kepO = kepler.from_statevector(statevector(r0,vi), body.GM, t0)
    kepI = kepler.from_statevector(statevector(r0,v0), body.GM, t0)
    taf = kepI.true_anomaly_by_distance(body.soi)
    return kepO, kepI, t0 - kepI.time_anomaly_by_ta(taf), t0


def solve_transfer(bodyA, rpe0, t0, bodyB, rpe3, t3, sun, refine=True, refine_maxiter=10):    
    stvi = bodyA.statevector(t0)
    stvf = bodyB.statevector(t3)
    KepT = kepler.lambert(stvi.r, t0, stvf.r, t3, sun.GM)
    
    def solve_eject_insert(kepT, t0, rpe0, bodyA, t3, rpe3, bodyB, sun):
        stv0 = kepT.statevector_by_time(t0)
        stv3 = kepT.statevector_by_time(t3)
        kepI, kepA, t0, t1 = solve_ejection_prograde(bodyA, stv0, t0, rpe0, sun)
        kepF, kepB, t2, t3 = solve_insertion_prograde(bodyB, stv3, t3, rpe3, sun)
        stv1 = kepT.statevector_by_time(t1)
        stv2 = kepT.statevector_by_time(t2)
        return (stv1, t1), (stv2, t2), (kepI, kepA, kepB, kepF)
    
    (stv1, t1), (stv2, t2), (kepI, kepA, kepB, kepF) = solve_eject_insert(kepT, t0, rpe0, bodyA, t3, rpe3, bodyB, sun)
        
    dv_eject = kepA.statevector_by_time(t0).v - kepI.statevector_by_time(t0).v
    dv_insert = kepF.statevector_by_time(t3).v - kepB.statevector_by_time(t3).v
    return [(None, None, kepI, bodyA)
            (t0, dv_eject, kepA, bodyA),
            (t1, None, kepT, sun),
            (t2, None, kepB, bodyB),
            (t3, dv_insert, kepF, bodyB)]

