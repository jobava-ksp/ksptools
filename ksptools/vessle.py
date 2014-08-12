from numpy import cross

from . import part
from . import body
from .util import unit, unitk, reject

class Stage(object):
    def __init__(self, parts, activeset, next):
        self.next = next                 # next stage
        self.parts = parts               # set of physically attached parts
        self.activeset = activeset       # set of active parts
    
    def total_mass(self):
        return self.part.mass + sum(n.total_mass() for n in self.next)
    
    def total_dragcoef(self):
        return self.part.dragcoef + sum(n.total_dragcoef() for n in self.next)
    
    def mass(self):
        raise sum(p.mass for p in self.parts)
    
    def dragcoef(self):
        return sum(p.dragceof() for p in parts)
    
    def thrust_func(self):
        def func(throttle, atm):
            return unitk * sum(e.thrust(throttle, atm) for e in self.activeset)
    
    def start(self):
        raise NotImplementedError


class Vessle(body.Body):
    def __init__(self, kename, name, state, stages):
        body.Body.__init__(self, keyname, name, state, orientation, mass=0.0)
        self.stages = stages
        self.orientation = orientation
        self.throttle = 0.5
        
    def prestep(self, world_time, dt):
        self.mass = stages[-1].total_mass()
        self.dragcoef = stages[-1].total_dragcoef(
        #...
    
    def poststep(self, world_time, dt):
        pass
    
    def bodyforce(self, world_time, dt):
        Ft = self.stages[-1].thrust_func()
        return self.orientation * Ft

