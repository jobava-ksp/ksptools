import copy

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


class ResourceTank(object):
    def __init__(self, resourcetype, quantity):
        self.resourcetype = resourcetype
        self.quantity = quantity
    
    @property
    def mass(self):
        return self.quantity * self.resourcetype.unitdensity
    
    def drain(self, quantity=None):
        if quantity is None:
            self.quantity = 0
        else:
            self.quantity += quantity
    
    def fill(self, quantity):
        self.quantity += quantity


class PartType(object):
    def __init__(self, name, partclass, title, cost, mass, coefdrag, resources):
        self.name = name
        self.partclass = partclass
        self.title = title
        self.cost = cost
        self.mass = mass
        self.coefdrag = coefdrag
        self.resources = dict((r.resourcetype.name, copy.copy(r)) for r in resources)


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
    
    def isp(self, atm=0.):
        return self.ispvac + self.ispscale*max(0, min(1, atm))
    
    def ff(self, scale=1., atm=0.):
        return (scale * self.maxthrust) / (9.80665 * self.isp(atm))
    
    def thrust(self, scale=1., atm=0.):
        return scale * self.maxtrhust


class FueledEnginePartType(EnginePartType):
    def __init__(self, name, title, cost, drymass, wetmass, coefdrag, fueltype, maxthrust, isp, ispatm, ff, ffatm):
        res = ResourceTank(fueltype, wetmass - drymass)
        PartType.__init__(self, name, ['sink', 'source'], title, cost, drymass, coefdrag, [res])
        self.resourcetype = fueltype
        self.maxthrust = maxthrust
        self.ispvac = isp
        self.ispatm = ispatm
        self.ffvac = ff
        self.ffatm = ffatm
        self.ispscale = ispatm - isp


class FuelTankPartType(PartType):
    def __init__(self, name, title, cost, drymass, wetmass, coefdrag, resourcetype):
        res = ResourceTank(resourcetype, wetmass - drymass)
        PartType.__init__(self, name, ['source'], title, drymass, coefdrag, [res])


class ChutePartType(PartType):
    def __init__(self, name, title, cost, mass, coefdrag, semidrag, fulldrag):
        PartType.__init__(name, ['chute'], title, cost, mass)
        self.coefdrag = coefdrag
        self.semidrag = semidrag
        self.maxdrag = maxdrag


class Part(object):
    def __init__(self, parttype, resource_quantity_mult=dict()):
        for k, v in vars(parttype):
            setattr(self, k, v)
        for k, f in resource_quantity_mult:
            self.resources[k].quantity = f * self.resources[k].quantity


class EnginePart(Part):
    def __init__(self, parttype):
        Part.__init__(self, parttype)
        self.isp = parttype.isp
        self.ff = parttype.ff
        self.thrust = parttype.thrust


class PartGroup(Part):
    def __init__(self, parts):
        self.parts = parts
        self.jointanks()
    
    def jointanks(self):
        self.resources = dict()
        for p in self.parts:
            for k, tank in p.resources:
                if k not in self.resources:
                    self.resources[k] = copy.copy(tank)
                else:
                    self.resources[k].quantity += tank.quantity
                tank.quantity = 0
    
    @property
    def mass(self):
        return sum(p.mass for p in self.parts) + sum(tank.mass for tank is self.resources.values())
    

def group(parts):
    if isinstance(parts, PartGroup):
        return parts
    else:
        return PartGroup(parts)


