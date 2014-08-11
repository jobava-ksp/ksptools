import re
import numpy as np
import numpy.linalg as npla

from . import orbit
from . import atmosphere

from .locallity import OrbitalState, FixedState
from .util import unit

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
        self._state = state
        self.sattelites = set()
    
    def __eq__(self, other):
        return self.keyname == other.keyname
    
    def __hash__(self):
        return hash(self.keyname)
    
    @property
    def kepler(self):
        return self._state.kepler
    
    @property
    def velocity(self):
        return self._state.velocity
    
    @property
    def position(self):
        return self._state.position
    
    @property
    def global_position(self):
        p = self._state.refbody
        r = self._state.position
        if p is not None:
            return r + p.global_position
        else:
            return r
    
    @property
    def global_velocity(self):
        p = self._state.refbody
        v = self._state.velocity
        if p is not None:
            return v
        else:
            return v + p.global_velocity
    
    @property
    def parent(self):
        return self._state.refbody
    
    @property
    def period(self):
        return self._state.period
    
    @property
    def surfnorm(self):
        return unit(self._state.position)
    
    def envforce(self, r, v, nbody=False, system=None):
        if nbody is None:
            Fg = (-unit(r))*(self.mass * self.refbody.std_g_param)/np.dot(r,r)
        else:
            raise NotImplementedError
        if self.refbody.atmposphere is not None and npla.norm(r) < (self.eq_radius + self.atmosphere.height):
            Fd = (-unit(v))*(0.5*dot(v,v)*self.coefDragArea()
        else:
            Fd = np.array([0.,0.,0.])
        return Fd + Fg
    
    def step(self, world_time, dt, nbody=False, system=None):
        r0, v0 = self._state.rv(self.epoch)
        envforce = self.envforce(r0, v0, nbody, system)
        bodyforce = self.bodyforce(world_time, dt)
        a = (envforce + bodyforce) / self.mass
        r = r0 + v0*dt + a*(dt**2)/2
        v = v0 + a*dt
        self._state = type(self._state).from_rv(self.refbody, r, v, self.body, self.epoch + dt)
    
    def coefDragArea(self):
        raise NotImplementedError
    
    def bodyforce(self, world_time, dt):
        raise NotImplementedError
    
    def prestep(self, world_time, dt):
        raise NotImplementedError
    
    def poststep(self, world_time, dt):
        raise NotImplementedError


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
            return OrbitalState(parent, orbit.KeplerOrbit.from_config(conf_parser, section, parent.std_g_param), None)
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
            if b.parent is None:
                self.centers.add(b)
            b.system = self
    
    def __getitem__(self, key):
        return self.bodies[key]
    
    def __setitem__(self, key, val):
        self.bodies[key] = val
        
