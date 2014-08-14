from numpy import cross

from . import part
from . import body
from .util import unit, unitk, reject

class AssemblyStage(part.PartGroup):
    def __init__(self, parts=[], name=None):
        part.PartGroup.__init__(self, parts, name)
        self.activation_stages = [part.PartGroup([]), part.PartGroup([])]
    
    def setstages(self, part_groups):
        self.activation_stages = [p if p is not None else part.PartGroup([]) for p in part_groups]
    
    def _pushignition(self, ignition_group):
        self.activation_stages[-1].addpart(ignition_group)
    
    def _popignition(self):
        ign = self.activation_stages[-1]
        self.activation_stages[-1] = part.PartGroup([])
        return ign
    
    def _popsep(self):
        sep = self.activation_stages[0]
        self.activation_setages = self.activation_stages[1:]
        return sep
    
    def autolinkengines(self):
        for engine in self.partsbyclass('engine'):
            if engine.tank is None:
                engine.link(self)


class Assembly(object):
    def __init__(self, assembly_stages):
        self.assembly_stages = assembly_stages
    
    @property
    def first(self):
        return self.assembly_stages[-1]
    
    @property
    def last(self):
        return self.assembly_stages[0]
    
    @staticmethod
    def stack(top, bottom):
        top.last._pushignition(bottom.first._popsep())
        return Assembly(top.assembly_stages + bottom.assembly_stages)
    
    @staticmethod
    def booster(core, stage, feed_part=None):
        stage.last._pushignition(core.last._popignition())
        if feed_part is not None:
            feed_part.link(stage, core)
            stage.addpart(feed_part)
            core.activation_stages[-1].addpart(feed_part)
        return Assembly.stack(core, stage)


def single_assembly(part_list, name=None, stages=['seperation', 'ignition']):
    stage = AssemblyStage(part_list, name)
    stage.setstages(stage[s] for s in stages)
    return Assembly([stage])

stack_assembly = Assembly.stack
booster_assembly = Assembly.booster


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

