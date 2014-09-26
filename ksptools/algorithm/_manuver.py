from . import _ode

from .._kepler import KeplerOrbit as kepler
from .._vector import statevector
from .._math import unitk, unit
from ..path._path import Path as path

from numpy import array, arccos, dot, cross, imag, pi, roots, sqrt
from numpy.linalg import norm

#         #
# Hohmann #
#         #

def solve_hohmann_transfer(body, kepA, ts, tpe, tap=None):
    ## Start at periapsis or apopasis
    stvA_pe, tA_pe = kepA.statevector_by_ta(0), kepA.time_at_ta(0, ts)
    stvA_ap, tA_ap = kepA.statevector_by_ta(pi), kepA.time_at_ta(pi, ts)    
    
    if tap is None or tpe == tap:
        ## possible transfers to a circular orbit
        th = sqrt(2*body.GM)*sqrt(tpe**2/(2*tpe))
        transfer_list = [
            (stvA_pe, tA_pe, tpe, tpe, body.GM),
            (stvA_ap, tA_ap, tpe, tpe, body.GM)]
    else:
        ## possible transfers to an elliptical orbit
        th = sqrt(2*body.GM)*sqrt((tap*tpe)/(tap + tpe))
        transfer_list = [
            (stvA_pe, tA_pe, tpe, tap, body.GM),
            (stvA_pe, tA_pe, tap, tpe, body.GM),
            (stvA_ap, tA_ap, tpe, tap, body.GM),
            (stvA_ap, tA_ap, tap, tpe, body.GM)]
    
    def _hohmann_dv(r0, v0, ta, tb, u):
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
    vTi = cross(kepA.normal, unit(stvA.r)*(hT/norm(stvA.r))
    stvTi = statevector(stvA.r, vTi)
    kepT = kepler.from_statevector(stvTi, u, tA)
    tB = tA + kepT.period/2.0
    stvTf = kepT.statevector_by_time(tB)
    
    ## final orbit (kepB)
    stvB = statevector(stvTf.r, unit(stvTf.v)*(hB/norm(stvTf.r)))
    kepB = kepler.from_statevector(stvB, u, tB)
    
    return path([
            path.coast(body, kepA, ts, tA),
            path.impulse(body, stvA, stvTi, tA),
            path.coast(body, kepT, tA, tB),
            path.impulse(body, stvTf, stvB, tB),
            path.coast(body, kepB, tB)])

#                        #
# Bi-elliptic Rondezvous #
#                        #

def solve_bielliptic_rondezvous(body, kepA, kepB, ts):
    rel_an = unit(cross(kepA.normal, kspB.normal)) # relative ascending node
    rel_dn = -unit(rel_dn)                         # relative descending node
    
    def solve_at_node(node, n, m):
        thetaA = kepA.true_anomaly_by_vector(node)
        tA = kepA.time_at_ta(thetaA, ts) + n * kepA.period
        rvA = kepA.statevector_by_ta(thetaA)
        rA = norm(rvA.r)
        
        thetaB = kepB.true_anomaly_by_vector(node)
        tB = kepB.time_at_ta(thetaB, tA) + m * kepB.period
        rvB = kepB.statevector_by_ta(thetaB)
        rB = norm(rvB.r)
        
        if tB >= tA:
            return []
        
        coefs = array([
            2,
            3*(rA + rB),
            3*(rA**2 + rB**2),
            rA**3 + rB**3 - 8*body.GM*(tB-tA)])
        return [((tA, rA, rvA, tB, rB, rvB, rT) for rT in roots(coefs) if rT > 0 and imag(rT) == 0]
    
    for n in itertools.count():
        for m in range(int(kepB.period/kepA.period) + 1):
            for tA, rA, rvA, tB, rB, rvB, rT in itertools.chain(solve_at_node(rel_an, rel_dn)):
                hT0 = sqrt(2*body.GM)*(rA*rT/(rA + rT))
                hT1 = sqrt(2*body.GM)*(rT*rB/(rT + rB))
                tT = (rA + rT)**3/(16*body.GM) + tA
                vT0i = cross(kepA.normal, rvA.r/rA)*(hT0/rA)
                vT1f = cross(kepB.normal, rvB.r/rB)*(hT1/rB)
                stvT0i = statevector(rvA.r, vT0i)
                stvT1f = statevector(rvB.r, vT1f)
                kepT0 = kepler.from_statevector(stvT0i, body.GM, tA)
                kepT1 = kepler.from_statevector(stvT1f, body.GM, tB)
                stvT0f = kepT0.statevector_by_time(tT)
                stvT1i = kepT1.statevector_by_time(tT)
                
                yield path([
                    path.coast(body, kepA, ts, tA),
                    path.impulse(body, rvA, stvT0i, tA),
                    path.coast(body, kepT0, tA, tT),
                    path.impulse(body, stvT0f, stvT1i, tT),
                    path.coast(body, kepT1, tT, tB),
                    path.impulse(body, stvT1f, rvB, tB),
                    path.coast(body, kepB, tB)])

#      #
# Burn #
#      #

def solve_burn(body, kepA, kepB, tc, isp, m0, dm):
    gafunc = _ode.gravity_accel_func(body.GM)
    dv = kepB.statevector_by_time(tc).v - kepA.statevector_by_time(tc).v
    tafunc, dmfunc = _ode.thrust_accel_func(unit(dv), isp, dm)
    
    def func(x):
        t0, delt = x
        stv0 = kepA.statevector_by_time(t0)
        stv1 = _ode.simrvm(stv0, m0, t0, t0+delt, lambda r,v,t,m: gafunc(r,v,t) + tafunc(r,v,t,m), dmfunc)
        return norm(kebB.statevector_by_time(t0 + delt)._vector - stv1._vector)
    
    ti, delt = minimize(func, [tc, 0]).x
    return path.burn(body, kepA, kepB, ti, ti+delt)


def solve_ejection_prograde(body, stv0, t0, rpe, sun=None):
    if sun is None:
        sun = body.parent_node
    
    ### orbital parameters (perifocal) ###
    vesc = sun.relative_statevector(body, stv0, t0).v
    spec_energy = dot(vesc,vesc)/2 - body.GM/body.soi
    sma = -body.GM/(2*spec_energy)
    e = 1 + rpe * dot(vesc, vesc)/body.GM
    p = sma * (1-e**2)
    vpe = sqrt(dot(vesc, vesc) + 2*body.GM/rpe)
    
    ### orientation ###
    ta = arccos(((p/rpe)-1)/e) # angle from apse line / true anomaly
    B = arccos(1/e)            # angle of assymptote
    normal = unit(corss(cross(vesc, unitk), vesc))
    r0 = rpe * dot(rotaxis(normal, B), (-unit(vesc))).A1
    v0 = vpe * unit(cross(noraml, r0))
    
    kepA = kepler.from_statevector(statevector(r0, v0), body.GM, t0)
    tf = kepA.true_anomaly_by_distance(


