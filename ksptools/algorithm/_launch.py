from numpy import arcsin, array, cos, dot, mat, pi, sin
from numpy.linalg import norm
from scipy.constants import g as g0

from .._math import unit
from .._kepler import KeplerOrbit as kepler
from ._ode import StagedController


class _Phase(object):
    def __init__(self, fa=None, twr=None):
        self.fa = fa
        self.twr = twr
    
    def objfunc(self, body, state):
        raise NotImplementedError


class _AltPhase(_Phase):
    def __init__(self, alt, fa=pi/2, twr=None):
        _Phase.__init__(self, fa, twr)
        self.tgt_alt = alt
    
    def objfunc(self, body, state):
        rmv, ijkllav, pav, t = state
        return (self.tgt_alt - ijkllav[5])**2


class _LaunchController(StagedController):
    def __init__(self):
        StagedController.__init__(self)
    
    def _staged_update(self, stage, depleted, maxtwr, state):
        mrv, ijkllav, pav, t = state
        m,r,_ = mrv
        i,j,k,_,_,_,_ = ijkllav
        _,_,_ = pav
        IJK = mat([i,j,k])
        
        ## get phase info... (fa, twr, body)
        fa, twr, body = self.tgt_fa, self.tgt_twr, self.body
        
        
        ## adjust twr ##
        if twr is not None:
            twr = min(1, twr)
        else:
            twr = maxtwr
        
        ## adjust fa such that Fg <= Tz
        if twr > 0:
            fa = max(0, max(arcsin(1/twr), fa))
        
        ## adjust throttle ##
        if twr > 0:
            f = min(1, twr/maxtwr)
        else:
            f = 0
        return f, dot(IJK.I, array([cos(fa), 0, sin(fa)])).A1, (depleted and twr > 0)
    
    def dophase(self, phase, stage0, t0, stv0, body, **kwargs):
        objfunc = phase.objfunc
        self.tgt_fa = phase.fa
        self.tgt_twr = phase.twr
        self.body = body
        return self.sim_to_objective(stage0, stv0, t0, body, objfunc, **kwargs)

def altphase(alt):
    return _AltPhase(alt)

def launch(phases, stage, t0, stv0, body):
    controller = _LaunchController()
    for phase in phases:
        stv0, stage, t0 = controller.dophase(phase, stage, stv0, t0, body)
    return stage, t0, stv0

