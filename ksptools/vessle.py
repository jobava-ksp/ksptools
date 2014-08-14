from numpy import cross

from . import part
from . import body
from .util import unit, unitk, reject

class Assembly(object):
    def __init__(self, prev_assembly=None):
        self.prev_assembly = prev_assembly
        self.engines = []
        
        self.part_group = part.PartGroup([])
        if prev_assembly is not None:
            self.assembly_group = part.group(None, self.part_group, prev_assembly.assembly_group)
            self.stages = prev_assembly.stages
        else:
            self.assembly_group = self.part_group
            self.stages = [None]
    
    def linkengines(self):
        for engine_group in self.engines:
            for engine in engine_group.partsbyclass('engine'):
                if engine.tank is None:
                    engine.link(self.part_group)
    
    def addstage(self, activationset):
        self.stages += [Stage(self.assembly_group, self.stages[-1], activationset)]
    
    def addpart(self, part):
        self.part_group.parts.append(part)
        return self
    
    def addengine(self, part):
        self.addpart(part)
        self.engines.append(part)
        return self
    
    def addchute(self, part):
        self.addpart(part)
        self.addstage([part])
        return self
    
    def addstack(self, seppart):
        self.linkengines()
        self.addstage([seppart] + self.engines)
        assembly = Assembly(self)
        assembly.addpart(seppart)
        return assembly
    
    def addbooster(self, seppart, feedpart):
        self.linkengines()
        self.addstage([seppart] + [feedpart] + self.engines)
        
        feedpart.link(assembly.part_group, self.part_group)
        assembly = Assembly(self)
        assembly.addpart(seppart)
        assembly.addpart(feedpart)
        assembly.engines += self.engines
        return assembly
    
    def finish(self):
        self.linkengines()
        self.addstage(self.engines)
        return self

class Stage(object):
    def __init__(self, partgroup, next, activationset):
        self.partgroup = partgroup
        self.next = next
        self.activationset = activationset
    
    @property
    def mass(self):
        return self.partgroup.mass
    
    @property
    def dragcoef(self):
        return self.partgroup.dragcoef
    
    def start(self):
        for part in self.activationset:
            part.activate()


class Vessle(body.Body):
    def __init__(self, kename, name, state, stages):
        body.Body.__init__(self, keyname, name, state, orientation, mass=0.0)
        self.stages = stages
        self.orientation = orientation
        self.throttle = 0.5
        
    def prestep(self, world_time, dt):
        pass
    
    def poststep(self, world_time, dt):
        pass
    
    def bodyforce(self, world_time, dt):
        pass

