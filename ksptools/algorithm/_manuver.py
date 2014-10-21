import itertools

from . import _ode
from .._kepler import KeplerOrbit as kepler
from .._vector import statevector
from .._math import rotaxis, unit, unitk
#from ..path._path import Path as path

from numpy import array, arccos, dot, cross, imag, pi, roots, sqrt
from numpy.linalg import norm
from scipy.optimize import minimize

#         #
# Hohmann #
#         #

def solve_hohmann_transfer(body, kepA, ts, tgt_pe, tgt_ap=None):
    """Solve a hohmann transfer from one orbit to another coaxil orbit.
    
    @param body The central body.
    @param KepA The starting orbit.
    @param ts Starting time for solutions. All solutions occur after ts.
    @param tgt_pe The periapsis of the target orbit, or the radius for a circular target orbit.
    @param tgt_ap The apoapsis of the target orbit, or None if circular. (default is None)
    """
    ## Starting points at periapsis and apopasis
    stvA_pe, tA_pe = kepA.statevector_by_ta(0), kepA.time_at_ta(0, ts)
    stvA_ap, tA_ap = kepA.statevector_by_ta(pi), kepA.time_at_ta(pi, ts)    
    
    ## list of possiple transfers
    if tgt_ap is None or tgt_pe == tgt_ap:
        ## possible transfers to a circular orbit
        th = sqrt(2*body.GM)*sqrt(tgt_pe**2/(2*tgt_pe))
        transfer_list = [
            (stvA_pe, tA_pe, tgt_pe, tgt_pe, body.GM),
            (stvA_ap, tA_ap, tgt_pe, tgt_pe, body.GM)]
    else:
        ## possible transfers to an elliptical orbit
        th = sqrt(2*body.GM)*sqrt((tgt_ap*tgt_pe)/(tgt_ap + tgt_pe))
        transfer_list = [
            (stvA_pe, tA_pe, tgt_pe, tgt_ap, body.GM),
            (stvA_pe, tA_pe, tgt_ap, tgt_pe, body.GM),
            (stvA_ap, tA_ap, tgt_pe, tgt_ap, body.GM),
            (stvA_ap, tA_ap, tgt_ap, tgt_pe, body.GM)]
    
    def _hohmann_dv(r0, v0, ta, tb, u):
        '''Calculate dv for a hohmann transfer from ta to tb'''
        th = sqrt(2*u)*sqrt((ta*tb)/(ta + tb))
        hh = sqrt(2*u)*sqrt((ta*r0)/(ta + r0))
        
        dv0 = abs(v0 - hh/r0)
        dv1 = abs(hh/ta - th/ta)
        return dv0 + dv1
    
    ## find the best transfer
    fdv = lambda x: _hohmann_dv(norm(x[0].r), norm(x[0].v), x[2], x[3], x[4])
    stvA, tA, r0, r1, u = min(transfer_list, key=fdv)
    
    ## orbital energies for transfer orbit (T) and final orbit (B)
    hT = sqrt(2*u)*sqrt((r0*norm(stvA.r))/(r0+norm(stvA.r)))
    hB = sqrt(2*u)*sqrt((r0*r1)/(r0+r1))
    
    ## Transfer orbit (kepT) and state vectors at tA and tB (stvTi, stvTf)
    vTi = cross(kepA.normal, unit(stvA.r)*(hT/norm(stvA.r)))
    stvTi = statevector(stvA.r, vTi)
    kepT = kepler.from_statevector(stvTi, u, tA)
    tB = tA + kepT.period/2.0
    stvTf = kepT.statevector_by_time(tB)
    
    ## final orbit (kepB)
    stvB = statevector(stvTf.r, unit(stvTf.v)*(hB/norm(stvTf.r)))
    kepB = kepler.from_statevector(stvB, u, tB)
    
    return [(tA, stvTi.v - stvA.v, kepT, body),
            (tB, stvB.v - stvTf.v, kebB, body)]


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

