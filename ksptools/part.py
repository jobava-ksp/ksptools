import itertools


class Catalogue(self):
    def add(self, parttype):
        setattr(self, parttype.name, parttype)
    
    def __getitem__(self, name):
        return getattr(self, name)


class ResourceType(object):
    def __init__(self, name, title, unitcost, unitdensity, flowmode):
        self.name = name
        self.title = title
        self.unitcost = unitcost
        self.unitdensity = unitdensity
        self.flowmode = flowmode
    
    def __call__(self, quantity):
        return ResourceTank(self, quantity)


class ResourceTank(object):
    def __init__(self, resourcetype, quantity):
        self.resourcetype = resourcetype
        self.quantity = quantity
    
    def copy(self):
        return ResourceTank(self.resourcetype, self.quantity)
    
    @property
    def mass(self):
        return self.quantity * self.resourcetype.unitdensity
    
    def drain(self, req_mass_amount):
        req_unit_amount = req_mass_amount / self.resourcetype.unitdensity
        unit_amount = min(self.quantity, req_unit_amount)
        self.quantity -= unit_amount
        return unit_amount * self.resourcetype.unitdensity, self.resourcetype, self.quantity > 0


class PartType(object):
    def __init__(self, name, partclass, title, radialsize, drycost, drymass, coefdrag, resources):
        self.name = name
        self.partclass = partclass
        self.title = title
        self.radialsize = radialsize
        self.drycost = drycost
        self.drymass = drymass
        self.coefdrag = coefdrag
        self.resources = dict((r.name, ResourceTank(r, q)) for q, r in resources)
    
    def __mul__(self, other):
        if isinstance(other, list):
            return map(self.__mul__, other)
        elif isinstance(other, tuple):
            quantity, name = other
            return self(q,n)
        return self(other)
    
    def __rmul__(self, other):
        return self.__mul__(other)
    
    def __call__(self, quantity=1, name=None):
        if quantity == 1:
            return type(self).instancetype()(self, name)
        else:
            return PartGroup([type(self).instancetype()(self, None) for i in range(quantity)], name)
    
    @property
    def mass(self):
        return self.drymass + sum(t.mass for t in self.resources.values())
    
    @property
    def coefdragArea(self):
        return self.mass * self.coefdrag * 0.008
    
    @classmethod
    def instancetype(cls):
        return Part


class Part(PartType):
    def __init__(self, parttype, name=None):
        for k, v in vars(parttype).items():
            if k not in ['resources']:
                setattr(self, k, v)
        if name is not None:
            self.name = name
        self.resources = dict((name, tank.copy()) for name, tank in parttype.resources.items())
        self.active = False
    
    @classmethod
    def instancetype(cls):
        return cls
    
    def drain(self, required_type):
        for required_amount, restype in required_type:
            yield self.resources[restype.name].drain(required_amount) 
    
    def activate(self):
        self.active = True
    
    def partsbyclass(self, partclass):
        if partclass in self.partclass:
            yield self
    

class PartGroup(Part):
    def __init__(self, parts, name=None):
        self.parts = parts
        self.name = name
        self.partclass = ['group']
        self.tanks = []
    
    def __call__(self, quantity=1, name=None):
        if quantity == 1:
            return PartGroup([p(1) for p in self.parts], name)
        else:
            return PartGroup([self(1) for i in range(quantity)], name)
    
    def link(self, tank):
        self.tanks.append(tank)
    
    def unlink(self, tank):
        self.links.remove(tank)
    
    def drain(self, required_type_tuple):
        remaining_amount, required_type = zip(*required_type_tuple)
        drained_amount = [0]*len(required_type)
        
        for part in itertools.chain(self.tanks, self.partsbyclass('tank')):
            remaining_amount, drained_amount = self._drainfrom(part, remaining_amount, drained_amount, required_type)
            if sum(remaining_amount) <= 0:
                return zip(drained_amount, required_type, [False]*len(required_type))
        return zip(drained_amount, required_type, [True]*len(required_type))
    
    def _drainfrom(self, part, remaining_amount, drained_amount, required_type):
        drained, restype, depleeted = zip(*(part.drain(zip(remaining_amount, required_type))))
        remaining_amount = tuple(total - delta for total, delta in zip(remaining_amount, drained))
        drained_amount = tuple(total + delta for total, delta in zip(drained_amount, drained))
        return remaining_amount, drained_amount
    
    def partsbyclass(self, partclass):
        for childpart in self.parts:
            for p in childpart.partsbyclass(partclass):
                yield p
    
    @property
    def mass(self):
        return sum(p.mass for p in self.parts)
    
    @property
    def coefdrag(self):
        return sum(p.coefdrag for p in self.parts)
    
    @property
    def drymass(self):
        return sum(p.drymass for p in self.parts)
    
    def activate(self):
        Part.activate(self)
        for p in self.parts:
            p.activate()
    
    def __getitem__(self, name):
        for p in self.parts:
            if 'group' in p.partclass:
                for subp in p[name]:
                    yield subp
            if p.name == name:
                yield p


def group(name=None, *parts):
    return PartGroup(list(parts), name)


class CrewPodPartType(PartType):
    def __init__(self, name, title, radialsize, drycost, drymass, ceofdrag, kerbals, torque, resources):
        PartType.__init__(self, name, ['pod', 'crew', 'sas'], title, radialsize, drycost, drymass, ceofdrag, resources)
        self.torque = torque
        self.maxkerbals = kerbals


class ProbePodPartType(PartType):
    def __init__(self, name, title, radialsize, drycost, drymass, ceofdrag, torque, resources):
        PartType.__init__(self, name, ['pod', 'crew', 'sas'], title, radialsize, drycost, drymass, ceofdrag, resources)
        self.torque = torque


class DecouplerPartType(PartType):
    def __init__(self, name, title, radialsize, drycost, drymass, coefdrag, ejection_momentum):
        PartType.__init__(self, name, ['decoupler'], title, radialsize, drycost, drymass, coefdrag, [])
        self.ejection_momentum = ejection_momentum


class FuelTankPartType(PartType):
    def __init__(self, name, title, radialsize, drycost, drymass, coefdrag, resources):
        PartType.__init__(self, name, ['tank'], title, radialsize, drycost, drymass, coefdrag, resources)


class EnginePartType(PartType):
    def __init__(self, name, title, radialsize, drycost, drymass, coefdrag, fueltypes, minthrust, maxthrust, isp, ispatm, resources=[]):
        PartType.__init__(self, name, ['engine'], title, radialsize, drycost, drymass, coefdrag, resources)
        self.fueltypes = fueltypes
        self.maxthrust = maxthrust
        self.minthrust = minthrust
        self.ispvac = isp
        self.ispatm = ispatm
        self.ispscale = ispatm - isp
        self.thrustscale = maxthrust - minthrust
        if len(self.resources):
            self.tank = self
        else:
            self.tank = None
        self.scale = 1.0
    
    @classmethod
    def instancetype(cls):
        return EnginePart


class EnginePart(Part):
    def isp(self, atm=0.):
        return self.ispvac + self.ispscale*max(0, min(1, atm))
    
    def ff(self, throttle, atm=0.):
        return self.thrust(throttle, atm) / (9.80665 * self.isp(atm))
    
    def thrust(self, trhottle, atm=0.):
        return self.scale *  (self.minthrust + throttle * self.thrustscale)
    
    def link(self, tankpart):
        self.tank = tankpart
    
    def step(self, throttle, dt, atm):
        if self.active:
            isp = self.ispvac + self.ispscale*max(0, min(1, atm))
            thrust = self.scale *  (self.minthrust + throttle * self.thrustscale)
            ff = thrust/(9.80665 * isp)
            fuel_required = [(mass_percent*ff*dt, fftype) for mass_percent, fftype in self.fueltypes]
            fuel_drained, _, depleeted = zip(*(self.tank.drain(fuel_required)))        
            thrust *= sum(fuel_drained)/ff
            return thrust, any(depleeted)
        else:
            return 0, False


class ChutePartType(PartType):
    def __init__(self, name, title, radialsize, drycost, drymass, coefdrag, semidrag, fulldrag, atm_semi, alt_full):
        PartType.__init__(self, name, ['chute'], title, radialsize, drycost, drymass, coefdrag, [])
        self.coefdrag = coefdrag
        self.atm_semideploy = atm_semi
        self.semidrag = semidrag
        self.alt_fulldeploy = alt_full
        self.fulldrag = fulldrag
    
    @classmethod
    def instancetype(cls):
        return ChutePart
        

class ChutePart(Part):
    def __init__(self, *args, **kwargs):
        Part.__init__(self, *args, **kwargs)
        self.chute_state = 'stowed'
        self.stowdrag = self.coefdrag
    
    def step(self, atm, alt_agl):
        if self.active:
            if self.chute_state == 'stowed' and atm >= self.atm_semideploy:
                self.coefdrag = self.semidrag
                self.chute_state = 'semi'
            if self.chute_state == 'semi' and alt_agl <= self.alt_fulldeply:
                self.coefdrag = self.fulldrag
                self.chute_state = 'full'
   
    def stow(self):
        self.active = False
        self.coefdrag = self.stowdrag
        self.chute_state = 'stowed'


class FeedPartType(PartType):
    def __init__(self, name, title, drycost, drymass, coefdrag):
        PartType.__init__(self, name, ['feed'], title, None, drycost, drymass, coefdrag, [])
    
    @classmethod
    def instancetype(cls):
        return FeedPart


class FeedPart(Part):
    def __init__(self, *args, **kwargs):
        Part.__init__(self, *args, **kwargs)
        self.frompart = None
        self.topart = None
    
    def link(self, frompart, topart):
        self.frompart = frompart
        self.topart = topart
        self.topart.link(self.frompart)
    
    def unlink(self):
        self.topart.unlink(self.frompart)
    
    def activate(self):
        Part.activate(self)
        self.unlink()

