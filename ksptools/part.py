import copy

class ResourceType(object):
    def __init__(self, name, title, unitcost, unitdensity, flowmode):
        self.name = name
        self.title = title
        self.unitcost = unitcost
        self.unitdensity = unitdensity
        self.flowmode = flowmode
    
    def __eq__(self, other):
        return self.name == other.name
    
    def __hash__(self):
        return hash(self.name)


class ResourceTank(object):
    def __init__(self, resourcetype, quantity):
        self.resourcetype = resourcetype
        self.quantity = quantity
    
    def copy(self):
        return ResourceTank(self.resourcetype, self.quantity)
    
    @property
    def mass(self):
        return self.quantity * self.resourcetype.unitdensity


class PartType(object):
    def __init__(self, name, partclass, title, radialsize, cost, mass, coefdrag, resources):
        self.name = name
        self.partclass = partclass
        self.title = title
        self.radialsize = radialsize
        self.cost = cost
        self.mass = mass
        self.coefdrag = coefdrag
        self.resources = dict((r.name, ResourceTank(r, q)) for q, r in resources)
    
    def new(self, quantity):
        return list(Part(self) for i in range(quantity))


class Part(object):
    def __init__(self, parttype):
        for k, v in vars(parttype):
            if k not in ['resources']:
                setattr(self, k, v)
        self.resources = dict((name, tank.copy()) for name, tank in parttype.resources)


class CrewPodPartType(PartType):
    def __init__(self, name, title, radialsize, cost, mass, ceofdrag, kerbals, torque, resources):
        PartType.__init__(self, name, ['pod', 'crew', 'sas'], title, radialsize, cost, mass, ceofdrag, resources)
        self.torque = torque
        self.maxkerbals = kerbals


class DecouplerPartType(PartType):
    def __init__(self, name, title, radialsize, cost, mass, coefdrag, ejection_momentum):
        PartType.__init__(self, name, ['decoupler'], title, radialsize, cost, mass, coefdrag, [])
        self.ejection_momentum = ejection_momentum


class FuelTankPartType(PartType):
    def __init__(self, name, title, radialsize, cost, mass, coefdrag, resources):
        PartType.__init__(self, name, ['tank'], title, radialsize, cost, mass, coefdrag, resources)


class EnginePartType(PartType):
    def __init__(self, name, title, radialsize, cost, mass, coefdrag, fueltypes, minthrust, maxthrust, isp, ispatm):
        PartType.__init__(self, name, ['engine'], title, radialsize, cost, mass, coefdrag, [])
        self.fueltypes = fueltypes
        self.maxthrust = maxthrust
        self.minthrust = minthrust
        self.ispvac = isp
        self.ispatm = ispatm
        self.ispscale = ispatm - isp
    
    def isp(self, atm=0.):
        return self.ispvac + self.ispscale*max(0, min(1, atm))
    
    def ff(self, throttle, scale=1., atm=0.):
        return scale * self.thrust(throttle, scale, atm) / (9.80665 * self.isp(atm))
    
    def thrust(self, trhottle, scale=1., atm=0.):
        return scale * (self.minthrust + throttle * (self.maxthrust - self.minthrust))


class ChutePartType(PartType):
    def __init__(self, name, title, radialsize, cost, mass, coefdrag, semidrag, fulldrag, atm_semi, alt_full):
        PartType.__init__(self, name, ['chute'], title, radialsize, cost, mass, coefdrag, [])
        self.coefdrag = coefdrag
        self.atm_semideploy = atm_semi
        self.semidrag = semidrag
        self.alt_fulldeploy = alt_full
        self.fulldrag = fulldrag


