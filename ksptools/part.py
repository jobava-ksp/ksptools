
class ResourceType(object):
    def __init__(self, name, title, unitcost, unitmass):
        self.name = name
        self.title = title
        self.unitcost = unitmass
        self.unitdensity = unitdensity
    
    def __eq__(self, other):
        return self.name == other.name
    
    def __hash__(self):
        return hash(self.name)


class ResourceSink(object):
    def __init__(self, resourcetype, quantity):
        self.resourcetype = resourcetype
        self.quantity = quantity
        self.sinks = []
        self.sources = []
    
    @staticmethod
    def combine(parts):
        res = dict()
        for part in parts:
            for r in part.resources:
                if r.resourcetype.name not in res:
                    res[r.resourcetype.name] = ResourceSink(r.resourcetype, 0)
                res[r.resource.name].quantity += r.quantity
        return res
    
    @staticmethod
    def link(sink, source):
        sink.sources.add(source)
        source.sinks.add(sink)


class PartType(object):
    def __init__(self, name, partclass, title, cost, mass, coefdrag, resources):
        self.name = name
        self.partclass = partclass
        self.title = title
        self.cost = cost
        self.mass = mass
        self.coefdrag = coefdrag
        self.resources = resources


class Part(object):
    def __init__(self, parttype):
        self.parttype = parttype
    
    def __eq__(self, other):
        


class EnginePartType(PartType):
    def __init__(self, name, title, cost, mass, coefdrag, fueltype, maxthrust, isp, ispatm, ff, ffatm):
        PartType.__init__(self, name, ['sink'], title, cost, mass, coefdrag, [])
        self.resourcetype = fueltype
        self.maxthrust = maxthrust
        self.ispvac = isp
        self.ispatm = ispatm
        self.ffvac = ff
        self.ffatm = ffatm
        self.ispscale = ispatm - isp
        self.sources = set()
    
    def isp(self, atm=0.):
        return self.ispvac + self.ispscale*max(0, min(1, atm))
    
    def ff(self, scale=1., atm=0.):
        return (scale * self.maxthrust) / (9.80665 * self.isp(atm))
    
    def trhust(self, scale=1., atm=0.):
        return scale * self.maxtrhust


class FueledEnginePartType(EnginePartType):
    def __init__(self, name, title, cost, drymass, wetmass, coefdrag, fueltype, maxthrust, isp, ispatm, ff, ffatm):
        res = ResourceSink(fueltype, wetmass - drymass)
        PartType.__init__(self, name, ['sink', 'source'], title, cost, drymass, coefdrag, [res])
        self.resourcetype = fueltype
        self.maxthrust = maxthrust
        self.ispvac = isp
        self.ispatm = ispatm
        self.ff = ff
        self.ffatm = ffatm
        self.ispscale = ispatm - isp
        self.sources = set()
        self.sinks = set()


class FuelTankPartType(PartType):
    def __init__(self, name, title, cost, drymass, wetmass, coefdrag, resourcetype):
        res = ResourceSink(resourcetype, wetmass - drymass)
        PartType.__init__(self, name, ['source'], title, drymass, coefdrag, [res])
        self.sinks = set()


class ChutePartType(PartType):
    def __init__(self, name, title, cost, mass, coefdrag, semidrag, fulldrag):
        PartType.__init__(name, ['chute'], title, cost, mass)
        self.coefdrag = coefdrag
        self.semidrag = semidrag
        self.maxdrag = maxdrag


