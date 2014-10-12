from numpy import arcsin, array, cos, dot, mat, pi, sin
from numpy.linalg import norm
from scipy.constants import g as g0

from .._math import unit
from .._kepler import KeplerOrbit as kepler
from ._ode import StagedController


class _Phase(object):
    def __init__(self, name, fa=None, twr=None):
        self.name = name
        self.fa = fa
        self.twr = twr
    
    def objfunc(self, body, state):
        raise NotImplementedError

    def oncomplete(self, stage, t, stv, body):
        lat, lon, alt, surf_v = body.llav(stv, t)
        print('t {:.1f}sec: {:.1f}m, {:.1f}m/s <{:.1f}m/s (horiz), {:+.1f}m/s (vert)>'.format(
            t, alt, norm(stv.v), norm(surf_v[0:2]), surf_v[2]))
        print('stage {}: {:.1f}kg, {:.1f}kg fuel, {:.1f}m/s stage dv, {:.1f}m/s total dv'.format(
            stage.name, stage.m0, stage.mp, stage.deltav(), stage.total_deltav()))
        print('{} deg N, {} deg E'.format(lat*(180/pi), lon*(180/pi)))
        print('--------')


class _AltPhase(_Phase):
    def __init__(self, alt, fa=pi/2, twr=None, name=''):
        _Phase.__init__(self, '{}: reach altitude {:.3e}m'.format(name, alt), fa, twr)
        self.tgt_alt = alt
    
    def objfunc(self, body, state):
        return (self.tgt_alt - state.alt)**2


class _ApoapsisPhase(_Phase):
    def __init__(self, ap, fa=pi/2, twr=None, name=''):
        _Phase.__init__(self, '{}: bring Apoapsis to {:.3e}m'.format(name, ap), fa, twr)
        self.tgt_ap = ap
    
    def objfunc(self, body, state):
        stv = state.stv
        _, _, _, ap = kepler.scalar_herpra(stv, body.GM)
        return (self.tgt_ap - ap)**2


class _PeriapsisPhase(_Phase):
    def __init__(self, pe, fa=0, twr=None, name=''):
        _Phase.__init__(self, '{}: bring Periapsis to {:.3e}m'.format(name, pe), fa, twr)
        self.tgt_pe = pe
    
    def objfunc(self, body, state):
        stv = state.stv
        _, _, pe, _ = kepler.scalar_herpra(stv, body.GM)
        return (self.tgt_pe - pe)**2


class _VPhase(_Phase):
    def __init__(self, v, fa=pi/2, twr=None, name=''):
        _Phase.__init__(self, '{}: reach {}m/s'.format(name, v), fa, twr)
        self.tgt_v = v
    
    def objfunc(self, body, state):
        return (self.tgt_v - norm(state.vsurf))**2


class _ReachApoapsisPhase(_Phase):
    def __init__(self, name=''):
        _Phase.__init__(self, '{}: reach apoapsis'.format(name), 0, 0)
    
    def objfunc(self, body, state):
        #_, _, _, ap = kepler.scalar_herpra(state.stv, body.GM)
        return -state.alt


class _FlightController(StagedController):
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
        #if twr > 0:
        #    fa = max(arcsin(1/twr), fa)
        
        ## adjust throttle ##
        if twr > 0:
            f = min(1, twr/maxtwr)
        else:
            f = 0
            
        return f, dot(IJK.I, array([cos(fa), 0, sin(fa)])).A1, (depleted and f > 0)
    
    def dophase(self, phase, stage0, t0, stv0, body, **kwargs):
        objfunc = phase.objfunc
        self.tgt_fa = phase.fa
        self.tgt_twr = phase.twr
        self.body = body
        return self.sim_to_objective(stage0, stv0, t0, body, objfunc, **kwargs)

def altphase(alt, fa=2/pi, twr=None, name=''):
    return _AltPhase(alt, fa, twr, name)

def apphase(ap, fa=pi/2, twr=None, name=''):
    return _ApoapsisPhase(ap, fa, twr, name)

def reachapphase(name=''):
    return _ReachApoapsisPhase(name)

def pephase(pe, fa=0, twr=None, name=''):
    return _PeriapsisPhase(pe, fa, twr, name)

def vphase(v, fa=pi/2, twr=None, name=''):
    return _VPhase(v, fa, twr, name)

def launch(phases, stage, t0, stv0, body):
    controller = _FlightController()
    for phase in phases:
        print(phase.name)
        stv0, stage, t0 = controller.dophase(phase, stage, stv0, t0, body)
        phase.oncomplete(stage, t0, stv0, body)
    return stage, t0, stv0

def stdlaunch(stage, t0, tgt_orbit=685.0e+3):
    from .._kerbolsystem import kerbolsystem
    body = kerbolsystem['kerbin']
    site = kerbolsystem['ksc']
    stv0 = site.statevector(t0)
    controller = _FlightController()
    phases = [
            altphase(10.0e+3),
            apphase(645.0e+3, pi/4),
            apphase(tgt_orbit, 0),
            reachapphase(),
            pephase(tgt_orbit)]
    stage, tf, stvf = launch(phases, stage, t0, stv0, body)
    return stage, tf, stvf, kepler.from_statevector(stvf, body.GM, tf)

