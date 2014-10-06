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
        return Stage(name, ms-(m0-m1), m0-m1, Tmax, isp_0, isp_1)


def _stage_link(stages):
    #stages = iter(stages)
    #top = next(stages)
    #for bottom in stages:
    #    bottom.next = top
    #    top = bottom
    #return bottom
    def linkfunc(a, n):
        n.next = a
        return n
    return functools.reduce(linkfunc, stages)


def _stage_copy(stage):
    return copy.deepcopy(stage)


class Stage(object):
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
        
    next = property(_get_next, _set_next)
    from_params = staticmethod(_stage_from_params)
    link = staticmethod(_stage_link)
    copy = staticmethod(_stage_copy)
    
    def deltav(self):
        return self.isp_0 * g0 * log(self.m0 / self.m1)
    
    def total_deltav(self):
        if self.next is not None:
            return self.deltav() + self.next.total_deltav()
        else:
            return self.deltav()
    
    def coefd(self, delta_m):
        if self.mp == 0:
            return self.coefd_ep
        else:
            return self.coefd_ep + delta_m*(self.coefd_e - self.coefd_ep)/self.mp
    
    def partial_deplete(self, mc):
        dm = self.m0 - mc
        return Stage(self.name, self.me, self.mp - dm, self.Tmax, self.isp_0, self.isp_1, self.coefd_e, self.coefd(dm), self.next)


class PartialStage(object):
    def __init__(self, me, mtfunc, Tmax, isp_0, isp_1, tankmass, tankstep, coefd_e=0.2, coefd_ep=0.3):
        self.me = me
        self.isp_0 = isp_0
        self.isp_1 = isp_1
        self.Tmax = Tmax
        self.coefd_e = coefd_e
        self.coefd_ep = coefd_ep
        self._tank_unit_mass = tankmass
        self._tank_unit_step = tankstep
        
    def build(self, name, mp):
        return Stage(name, self.me + self.mtfunc(mp), mp, self.Tmax, self.isp_0, self.isp_1, self.ceofd_e, self.ceofd_ep)
    
    def mtfunc(self, mp):
        units = int(self._tank_unit_step/mp) + 1
        return units * self._tank_unit_mass
    
    def gtfunc(self, mp):
        return self._tank_unit_mass


def minfuel(pl, names, partial_stages, dv, mintwr):
    from scipy.optimize import minimize
    plstage = Stage('pl', pl)
    
    def makestages(mp):
        stages = itertools.starmap(PartialStage.build, zip(partial_stages, names, mp))
        return Stage.link([plstage] + list(stages))
        
    def dvfunc(mp):
        stage0 = makestages(mp)
        return stage0.total_deltav() - dv
    
    def grad_dvfunc(mp):
        stage = makestages(mp)
        gmp = zeros(len(mp))
        for i in range(len(mp)):
            gmp[i] = stage.isp_0 * g0 * (partial_stages[i].gtfunc(mp[i])+1)/stage.m1
            stage = stage.next
        return gmp
    
    def fuelfunc(mp):
        return sum(mp)
    
    def grad_fuelfunc(mp):
        return array([1]*len(mp))
    
    mp0 = zeros(len(partial_stages))
    stage0 = plstage
    for i in reversed(range(len(partial_stages))):
        m0 = stage0.m0 + partial_stages[i].me
        mp0[i] = m0*(e**((dv/len(partial_stages))/(partial_stages[i].isp_0 * g0)) - 1)
        stage0 = Stage.link([stage0, PartialStage.build(names[i], mp0[i])])
    cos = ({'type': 'eq', 'func': dfvunc, 'jac': grad_dvfunc})
    minimize(fuelfunc, mp0, jac=grad_fuelfunc, constraints=cons)

