from numpy import cross

from . import part
from . import body
from .util import unit, unitk, reject

class Stage(object):
    def __init__(self, parts, next=[]):
        self.next = next
        self.part = part.group(parts)
    
    def total_mass(self):
        return self.part.mass + sum(n.total_mass() for n in self.next)
    
    def total_dragcoef(self):
        return self.part.dragcoef + sum(n.total_dragcoef() for n in self.next)
    
    def mass(self):
        raise NotImplementedError
    
    def dragcoef(self):
        raise NotImplementedError

class Vessle(body.Body):
    def __init__(self, kename, name, state, stages):
        body.Body.__init__(self, keyname, name, state, orientation, mass=0.0)
        self.stages = stages
        
    def prestep(self, world_time, dt):
        self.mass = stages[0].total_mass()
    
    def poststep(self, world_time, dt):
        pass
    
    def bodyforce(self, world_time, dt):
        

