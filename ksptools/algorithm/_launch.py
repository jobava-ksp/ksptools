from numpy import dot, mat, pi
from numpy.linalg import norm
from scipy.optimize import newton

from .._math import unit
from .._kepler import KeplerOrbit as kepler
from ._ode import StagedController


def target_rAp_func(target_rAp):
    def func(alt, rAp, tAp):
        return target_rAp - rAp

def target_alt_func(target_alt):
    def func(alt, rAp, tAp):
        return target_alt - alt


class LaunchPhaseController(StagedController):
    def __init__(self, tgt_dir, tgt_twr=None):
        StagedController.__init__(self)
        self.tgt_dir = tgt_dir
        self.tgt_twr = tgt_twr
    
    def _staged_update(self, stage_name, is_depleted, maxtwr, state):
        m_rv, body_ijk_llav, pav, t = state
        body, i, j, k, _, _, v_surf = body_ijk_llav
        
        IJK = mat([i,j,k])
        
        if self.tgt_twr is not None and maxtwr > 0:
            f = min(1, self.tgt_twr/maxtwr)
        else:
            f = 1
        
        return f, 0.2, dot(IJK, self.tgt_dir).A1, is_depleted
    
    def _solve(self, stages, stv0, body, t0, est_dt, objfunc):
        def func(dt):
            stv1, _, _ = self.sim(pl, stages, stv0, body, t0, dt)
            ## get altitude ##
            _, _, alt, _ = body.llav(stv1, t0+dt)
            ## get Ap ##
            kep = kepler.from_statevector(stv1, body.GM, t0+dt)
            rAp = norm(kep.statevector_by_ta(pi)
            tAp = kep.time_by_ta(pi, t0+dt)
            return objfunc(alt, rAp, tAp)
        
        dt = newton(func, est_dt)
        res = self.sim(pl, stages, stv0, body, t0, dt)
        return dt, res

