from numpy import log
from scipy.constants import g as g0


def _stage_from_params(name, ms, isp_0=0, isp_1=0, Tmax=0, Twr_0=0, Twr_1=0):
    if Twr_0 == 0:
        return Stage(name, ms)
    else:
        m1 = Tmax/(g0*Twr_1)
        m0 = Tmax/(g0*Twr_0)
        return Stage(name, ms-(m0-m1), m0-m1, Tmax, isp_0, isp_1)


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
    
    def deltav(self):
        return self.isp_0 * g0 * log(self.m0 / self.m1)
    
    def total_deltav(self):
        return self.isp_0 * g0 * log(self.m0 / self.m1)
    
    next = property(_get_next, _set_next)
    from_parameters = staticmethod(_stage_from_params)
    
    def coefd(self, delta_m):
        if self.mp == 0:
            return self.coefd_ep
        else:
            return self.coefd_ep + delta_m*(self.coefd_e - self.coefd_ep)/self.mp
    
    def partial_deplete(self, mc):
        dm = self.m0 - mc
        return Stage(self.name, self.me, self.mp - dm, self.Tmax, self.isp_0, self.isp_1, self.coefd_e, self.coefd(dm), self.next)



