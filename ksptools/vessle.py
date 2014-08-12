from numpy import cross

from . import part
from . import body
from .util import unit, unitk, reject

class Stage(object):
    def __init__(self, parts, engines=[], chutes=[], next=[]):
        self.next = next
        self.part = part.group(parts)
        self.engines = engines
        self.chutes = chutes
    
    def total_mass(self):
        return self.part.mass + sum(n.total_mass() for n in self.next)
    
    def total_dragcoef(self):
        return self.part.dragcoef + sum(n.total_dragcoef() for n in self.next)
    
    def mass(self):
        raise sum(p.mass for p in parts)
    
    def dragcoef(self):
        return sum(p.dragceof() for p in parts)
    
    def start(self):
        for e in self.engines:
            e.activate()
        for c in self.chutes:
            c.deploy()


class Vessle(body.Body):
    def __init__(self, kename, name, state, stages):
        body.Body.__init__(self, keyname, name, state, orientation, mass=0.0)
        self.stages = stages
        self.orientation = orientation
        
    def prestep(self, world_time, dt):
        self.mass = stages[-1].total_mass()
        #...
    
    def poststep(self, world_time, dt):
        raise NotImplementedError
    
    def bodyforce(self, world_time, dt):
    
    def dragcoef(self):
        return sum(s.dragcoef() for s in self.stages)

