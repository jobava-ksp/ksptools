from ..path._path import Path as path
from numpy import sqrt, dot, pi
from numpy.linalg import norm

#         #
# Hohmann #
#         #

def _hohmann_dv(r0, v0, ta, tb, u):
    th = sqrt(2*u)*sqrt((ta*tb)/(ta + tb))
    hh = sqrt(2*u)*sqrt((ta*r0)/(ta + r0))
    
    dv0 = abs(v0 - hh/r0)
    dv1 = abs(hh/ta - th/ta)
    return dv0 + dv1


def solve_hohmann_transfer(body, kepA, ts, tpe, tap=None):
    ## Start at periapsis or apopasis
    stvA_pe, tA_pe = kepA.statevector_by_ta(0), kepA.time_by_ta(0, ts)
    stvA_ap, tA_ap = kepA.statevector_by_ta(pi), kepA.time_by_ta(pi, ts)    
    
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
    
    ## find the best transfer
    fdv = lambda x: _hohmann_dv(norm(x[0].r), norm(x[0].v), x[2], x[3], x[4])
    stvA, tA, r0, r1, u = min(transfer_list, key=fdv)
    
    ## orbital energies for transfer orbit (T) and final orbit (B)
    hT = sqrt(2*u)*sqrt((r0*norm(stvA.r))/(r0+norm(stvA.r)))
    hB = sqrt(2*u)*sqrt((r0*r1)/(r0+r1))
    
    ## types ##
    vector_type = type(stvA)
    orbit_type = type(kepA)
    
    ## Transfer orbit (kepT) and state vectors at tA and tB (stvTi, stvTf)
    stvTi = vector_type(stvA.r, (stvA.v/norm(stvA.v))*(hT/norm(stvA.r)))
    kepT = orbit_type.from_statevector(stvTi, u, tA)
    tB = tA + kepT.period/2.0
    stvTf = kepT.statevector_by_time(tB)
    
    ## final orbit (kepB)
    stvB = vector_type(stvTf.r, (stvTf.v/norm(stvTf.v))*(hB/norm(stvTf.r)))
    kepB = orbit_type.from_statevector(stvB, u, tB)
    
    return path([
            path.coast(body, kepA, ts, tA),
            path.impulse(body, stvA, stvTi, tA),
            path.coast(body, kepT, tA, tB),
            path.impulse(body, stvTf, stvB, tB),
            path.coast(body, kepB, tB)])

#                        #
# Bi-elliptic Rondezvous #
#                        #

def solve_bielliptic_rondezvous(kepA, kepB):
    hA = kepA.specific_angular_momentum
    hB = kepB.specific_angular_momentum
    
