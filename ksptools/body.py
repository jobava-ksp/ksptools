import re

from . import orbit
from . import atmosphere


class State(object):
    def __init__(self, refbody, body=None, epoch=0.):
        self.refbody = refbody
        self.body = body
        self.epoch = epoch
    
    def asorbit(self):
        raise NotImplementedError
    
    def asvectors(self):
        raise NotImplementedError
    
    def rv(self, t):
        raise NotImplementedError


class OrbitalState(State):
    def __init__(self, refbody, kepler, body=None, epoch=0.):
        State.__init__(self, refbody, body, epoch)
        self.kepler = kepler
        self.rv = self.kepler.rv
        self.period = self.kepler.period
    
    def asorbit(self):
        return self
    
    def asvectors(self):
        r, v = self.kepler.rv(self.epoch)
        return VectorState(self.refbody, r, v, self.body, self.epoch)


class VectorState(State):
    def __init__(self, refbody, r, v, body=None, epoch=0.):
        State.__init__(self, refbody, body, epoch)
        self.velocity = v
        self.position = r
    
    def asorbit(self):
        kep = orbit.KeplerOrbit.from_rvu(
                self.velocity,
                self.position,
                self.refbody.std_g_param)
        return OrbitalState(self.refbody, kep, self.body, self.epoch)

    def asvectors(self):
        return self
    
    def rv(self, t):
        from numpy import isclose
        dt = t - epoch
        if all(isclose([dt],[0.0])):
            return self.position, self.velocity


class FixedState(State):
    def __init__(self, refbody=None, body=None, epoch=0.):
        State.__init__(self, refbody, body, epoch)


class Body(object):
    def __init__(self, keyname, name, state, mass=None, u=None):
        self.keyname = keyname
        self.name = name
        if u is None:
            self.std_g_param = mass * 6.67384e-11
            self.mass = mass
        else:
            self.std_g_param = u
            self.mass = u / 6.67384e-11
        self.state = state
        self.sattelites = set()
    
    def __eq__(self, other):
        return self.keyname == other.keyname
    
    def __hash__(self):
        return hash(self.keyname)
    
    def getorbit(self):
        return self.state.asorbit()
    
    def getvectors(self):
        return self.state.asvectors()


class CelestialBody(Body):
    def __init__(self, keyname, name, eq_radius, u, sidereal_rate, soi, body_atmosphere, body_orbit):
        Body.__init__(self, keyname, name, body_orbit, u=u)
        self.eq_radius = eq_radius
        self.sidereal_rate = sidereal_rate
        self.soi = soi
        self.atmosphere = body_atmosphere
    
    @classmethod
    def from_config(cls, conf_parser, section, system_dict):
        key_name = re.match(r'cbody\.(?P<name>[a-zA-Z0-9_]+)', section).group('name')
        orbit_section_name = 'cbody.{}.orbit'.format(key_name)
        atm_section_name = 'cbody.{}.atmosphere'.format(key_name)
        
        if conf_parser.has_option(section, 'parent'):
            parent_key = conf_parser.get(section, 'parent')
            parent = system_dict[parent_key]
        else:
            parent = None
        
        body_atmosphere = cls._atmosphere_from_config(conf_parser, atm_section_name)
        state = cls._orbit_from_config(conf_parser, orbit_section_name, parent)
        
        name = conf_parser.get(section, 'name')
        eq_radius = conf_parser.getfloat(section, 'radius')
        u = conf_parser.getfloat(section, 'gravitational_param')
        sidereal_rate = conf_parser.getfloat(section, 'sidereal_rate')
        soi = conf_parser.getfloat(section, 'soi')
        
        new_body = cls(key_name, name, eq_radius, u, sidereal_rate, soi, body_atmosphere, state)
        state.body = new_body
        if parent is not None:
            parent.sattelites.add(new_body)
        return new_body
    
    @classmethod
    def _orbit_from_config(cls, conf_parser, section, parent):
        if conf_parser.has_section(section):
            return OrbitalState(None, orbit.KeplerOrbit.from_config(conf_parser, section, parent.std_g_param), parent)
        else:
            return FixedState()
    
    @classmethod
    def _atmosphere_from_config(cls, conf_parser, section):
        if conf_parser.has_section(section):
            return atmosphere.Atmosphere.from_config(conf_parser, section)
        else:
            return None

class System(object):
    def __init__(self, bodies):
        self.bodies = dict((b.keyname, b) for b in bodies)
        self.centers = set()
        for b in bodies:
            if b.state.refbody is None:
                self.centers.add(b)
            b.system = self
    
    def __getitem__(self, key):
        return self.bodies[key]
    
    def __setitem__(self, key, val):
        self.bodies[key] = val
        
