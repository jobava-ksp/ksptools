import copy
import functools
import itertools
from numpy import array, e, log, zeros
from scipy.constants import g as g0


def _stage_from_params(name, ms, isp_0=0, isp_1=0, Tmax=0, Twr_0=0, Twr_1=0):
    if Twr_0 == 0:
        return Stage(name, ms)
    else:
        m1 = Tmax/(g0*Twr_1)
        m0 = Tmax/(g0*Twr_0)
        return _Stage(name, ms-(m0-m1), m0-m1, Tmax, isp_0, isp_1)


def _stage_link(stages):
    def linkfunc(a, n):
        n = _Stage.copy(n)
        n.next = a
        return n
    return functools.reduce(linkfunc, stages)


def _stage_expand(stage, unlink=True):
    def iterstages(s):
        while s is not None:
            if unlink:
                n = _Stage.copy(s)
                n.next = None
                yield n
            else:
                yield s
            s = s.next
    return list(reversed(list(iterstages(stage))))


def _stage_copy(stage):
    return copy.deepcopy(stage)


class _Stage(object):
    def __init__(self, name, me, mp=0, Tmax=0, isp_0atm=0, isp_1atm=0, coefd_e=0.2, coefd_ep=0.2, next=None):
        self.name = name
        self.me = me
        self.mp = mp
        self.Tmax = Tmax
        self.isp_0 = isp_0atm
        self.isp_1 = isp_1atm
        self.coefd_e = coefd_e
        self.coefd_ep = coefd_ep
        self.next = next
    
    def statevars(self):
        return self.name, self.m1, self.m0, self.Tmax, self.isp_0, self.isp_1, self.next
    
    def _get_next(self):
        return self._next
    
    def _set_next(self, next):
        self._next = next
        if next is not None:
            self.m1 = self.me + next.m0
            self.m0 = self.me + next.m0 + self.mp
        else:
            self.m1 = self.me
            self.m0 = self.me + self.mp
    
    def __rmul__(self, n):
        return _Stage(
                '{}x{}'.format(n, self.name),
                n*self.me,
                n*self.mp,
                n*self.Tmax,
                self.isp_0,
                self.isp_1,
                self.coefd_e,
                self.coefd_ep)
    
    def __add__(self, other):
        return _Stage.link([self, other])
    
    next = property(_get_next, lambda obj,val: type(obj)._set_next(obj, val))
    from_params = staticmethod(_stage_from_params)
    link = staticmethod(_stage_link)
    expand = staticmethod(_stage_expand)
    copy = staticmethod(_stage_copy)
    
    def deltav(self):
        return self.isp_0 * g0 * log(self.m0 / self.m1)
    
    def total_deltav(self):
        if self.next is not None:
            return self.deltav() + self.next.total_deltav()
        else:
            return self.deltav()
    
    def twr(self, g=g0):
        return self.Tmax/(self.m0*g)
    
    def coefd(self, delta_m):
        if self.mp == 0:
            return self.coefd_ep
        else:
            return self.coefd_ep + delta_m*(self.coefd_e - self.coefd_ep)/self.mp
    
    def partial_deplete(self, mc):
        dm = self.m0 - mc
        return _Stage(self.name, self.me, self.mp - dm, self.Tmax, self.isp_0, self.isp_1, self.coefd_e, self.coefd(dm), self.next)


class _BoosterStage(_Stage):
    def __init__(self, *args, **kwargs):
        _Stage.__init__(self, *args, **kwargs)
        self._init_args = (args, kwargs)
    
    def _set_next(self, next):
        if hasattr(self, '_next'):
            if next is not None:
                self_ff = self.Tmax/(g0*self.isp_0)
                self_dt = self.mp / self_ff
                self_dm = self.mp
                
                next_ff = next.Tmax/(g0*next.isp_0)
                next_dt = next.mp / next_ff
                next_dm = next.mp
                
                if self_dt < next_dt:
                    dm = min(next.mp, self_dt * next_ff)
                    next.mp -= dm
                    next.m0 -= dm
                    self.mp += dm
                    self.m0 = self.me + self.mp + next.m0
                    self.m1 = self.me + next.m0
                    self.isp_0 = sum(s.Tmax for s in [self, next])/sum(s.Tmax/s.isp_0 for s in [self, next])
                    self.isp_1 = sum(s.Tmax for s in [self, next])/sum(s.Tmax/s.isp_1 for s in [self, next])
                    self.Tmax += next.Tmax
                else:
                    
                    next_me = self.me + next.me
                    next_mp = self.mp - self_dm
                    next_m0 = next.m1 - next.me + (next_me + next_mp)
                    next_m1 = next.m1 - next.me + (next_me)
                    next_isp_0 = self.isp_0
                    next_isp_1 = self.isp_1
                    next_Tmax = self.Tmax
                  
                    self_me = 0
                    self_mp = next_mp + self_dm
                    self_m0 = next_m0 + self_me + self_mp
                    self_m1 = next_m0 + self_me
                    self_isp_0 = sum(s.Tmax for s in [self, next])/sum(s.Tmax/s.isp_0 for s in [self, next])
                    self_isp_1 = sum(s.Tmax for s in [self, next])/sum(s.Tmax/s.isp_1 for s in [self, next])
                    self_Tmax = self.Tmax + next.Tmax
                    
                    next.me, next.mp, next.m0, next.m1 = next_me, next_mp, next_m0, next_m1
                    next.isp_0, next.isp_1, next.Tmax = next_isp_0, next_isp_1, next_Tmax
                    self.me, self.mp, self.m0, self.m1 = self_me, self_mp, self_m0, self_m1
                    self.isp_0, self.isp_1, self.Tmax = self_isp_0, self_isp_1, self_Tmax
                    #print('bad')
                    '''
                    next.mp = 0
                    next.m0 = next.m1
                    self.mp += next.mp
                    self.m0 = self.me + self.mp + next.m0
                    self.m1 = self.me + next.m0
                    self.isp_0 = sum(s.Tmax for s in [self, next])/sum(s.Tmax/s.isp_0 for s in [self, next])
                    self.isp_1 = sum(s.Tmax for s in [self, next])/sum(s.Tmax/s.isp_1 for s in [self, next])
                    self.Tmax += next.Tmax
                    '''
                    
            else:
                raise NotImplementedError
        else:
            self.m0 = self.me
            self.m1 = self.me + self.mp
        self._next = next


class _MergedBoosterStage(_Stage):
    def __init__(self, *args, **kwargs):
        _Stage.__init__(self, *args, **kwargs)
        self._init_args = (args, kwargs)
    
    def _set_next(self, next):
        if hasattr(self, '_next'):
            if next is not None:
                self.isp_0 = sum(s.Tmax for s in [self, next])/sum(s.Tmax/s.isp_0 for s in [self, next])
                self.isp_1 = sum(s.Tmax for s in [self, next])/sum(s.Tmax/s.isp_1 for s in [self, next])
                self.Tmax += next.Tmax
            else:
                raise NotImplementedError
        else:
            self.m0 = self.me
            self.m1 = self.me + self.mp
        self._next = next


class _PartialStage(object):
    def __init__(self, name, me, Tmax, isp_0, isp_1, tmratio=0.125, tankstep=0.0625e+3, coefd_e=0.2, coefd_ep=0.2):
        self.name = name
        self.me = me
        self.isp_0 = isp_0
        self.isp_1 = isp_1
        self.Tmax = Tmax
        self.coefd_e = coefd_e
        self.coefd_ep = coefd_ep
        self._tank_mass_ratio = tmratio  # tank mass per fuel mass
        self._tank_unit_mass = tankstep
    
    def __rmul__(self, n):
        return type(self)(
                '{}x{}'.format(n, self.name),
                n*self.me,
                n*self.Tmax,
                self.isp_0,
                self.isp_1,
                self._tank_mass_ratio,
                n*self._tank_unit_mass,
                self.coefd_e,
                self.coefd_ep)
    
    def build(self, mp):
        return type(self)._stagetype(self.name, self.me + self.mtfunc(mp), mp, self.Tmax, self.isp_0, self.isp_1, self.coefd_e, self.coefd_ep)
    
    def burntime(self, mp):
        ff = self.Tmax/(self.isp_0*g0)
        return mp/ff
    
    def ff(self):
        return self.Tmax/(self.isp_0*g0)
    
    def mtfunc(self, mp):
        if self._tank_mass_ratio > 0:
            tank_mass = mp * self._tank_mass_ratio
            units = int(tank_mass / self._tank_unit_mass) + 1
            return units * self._tank_unit_mass
        else:
            return 0
    
    _stagetype = _Stage


class _PartialBoosterStage(_PartialStage):
    def __init__(self, *args, **kwargs):
        _PartialStage.__init__(self, *args, **kwargs)
    
    _stagetype = _BoosterStage


class _PartialMergedBoosterStage(_PartialStage):
    def __init__(self, *args, **kwargs):
        _PartialStage.__init__(self, *args, **kwargs)
    
    _stagetype = _MergedBoosterStage


def minimizefuel(pl, partial_stages, mintwr=None, dv=5000, fixed=[]):
    from scipy.optimize import minimize
    if mintwr is None:
        mintwr = [0]*len(partial_stages)
    plstage = _Stage('pl', pl)
    
    def makestages(mp, i=None):
        if i is None:
            i = len(partial_stages)-1
        stages = [_PartialStage.build(*t) for t in zip(partial_stages[:i+1], mp[:i+1])]
        return _Stage.link([plstage] + stages)
        
    def dvfunc(mp):
        stage0 = makestages(mp)
        sdv = stage0.total_deltav()
        return sdv - dv
    
    def twrfunc(i):
        def func(mp):
            stage = makestages(mp, i)
            tw = (stage.Tmax/(stage.m0*g0))
            return tw - mintwr[i]
        return func
    
    def fuelfunc(mp):
        return sum(mp)
    
    def grad_fuelfunc(mp):
        return array([1]*len(mp))
    
    ## initial value ##
    mp0 = zeros(len(partial_stages))
    stage0 = plstage
    for i in range(len(partial_stages)):
        m0 = stage0.m0 + partial_stages[i].me
        mp0[i] = m0*(e**((dv/len(partial_stages))/(partial_stages[i].isp_0 * g0)) - 1)
        stage0 = _Stage.link([stage0, _PartialStage.build(partial_stages[i], mp0[i])])
    
    ## constraints and bounds ##
    cons = []
    cons += [{'type': 'eq', 'fun': dvfunc, 'jac': None}]
    cons += [{'type': 'ineq', 'fun': twrfunc(i), 'jac': None} for i in range(len(mp0))]
    cons += [{'type': 'ineq', 'fun': lambda mp: mp[i] - minp, 'jac':None} for i, minp, _ in fixed if minp is not None]
    cons += [{'type': 'ineq', 'fun': lambda mp: maxp - mp[i], 'jac':None} for i, _, maxp in fixed if maxp is not None]
    bounds = [(0, None)]*len(mp0)
    
    ## fixed bounds ##
    for i, minp, maxp in fixed:
        bounds[i] = (minp, maxp)
        if minp is not None:
            mp0[i] = max(minp, mp0[i])
        if maxp is not None:
            mp0[i] = min(maxp, mp0[i])
    
    minres = minimize(fuelfunc, mp0, jac=grad_fuelfunc, bounds=bounds, constraints=cons, method='SLSQP', options={'maxiter':200})
    #return [(s.name, s.mp, s.twr(), s.deltav(), s.total_deltav()) for s in _Stage.expand(makestages(minres.x), False)[1:]]
    stage0 = plstage
    for i in range(len(minres.x)):
        stage = makestages(minres.x[:(i+1)])
        name = stage.name
        mp = minres.x[i]
        twr = stage.twr()
        dv = stage.deltav()
        sumdv = stage.total_deltav()
        yield name, mp, twr, dv, sumdv
    


