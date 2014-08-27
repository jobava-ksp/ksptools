import re
import numpy as np
import numpy.linalg as npla

from . import orbit
from . import atmosphere

from .locallity import OrbitalState, FixedState, GeocentricState
from .util import unit


class Body(object):
    def __init__(self, keyname, name, position, mass=None, u=None):
        """
        :type keyname: object
        :type name: str
        :type position: ksptools.locallity.State
        """
        self.keyname = keyname
        self.name = name
        if u is None:
            self.std_g_param = mass * 6.67384e-11
            self.mass = mass
        else:
            self.std_g_param = u
            self.mass = u / 6.67384e-11
        self.state = position
        self.sattelites = set()
    
    def __eq__(self, other):
        return self.keyname == other.keyname
    
    def __hash__(self):
        return hash(self.keyname)
    
    @property
    def parent(self):
        return self.state.refbody


class RigidBody(Body):
    def __init__(self, keyname, name, position, orientation, mass=0):
        Body.__init__(self, keyname, name, position, mass=mass)
        self.orientation = orientation


class CelestialBody(Body):
    def __init__(self, keyname, name, eq_radius, u, sidereal_rate, soi, body_atmosphere, body_orbit):
        """
        :type eq_radius: float
        :type u: float
        :type sidereal_rate: float
        :type soi: float
        :type body_atmosphere: ksptools.atmosphere.Atmosphere
        """
        Body.__init__(self, keyname, name, body_orbit, u=u)
        self.eq_radius = eq_radius
        self.sidereal_rate = sidereal_rate
        self.soi = soi
        self.atmosphere = body_atmosphere
    
    def geocentric(self, lon, lat, epoch, velocity=np.array([0,0,0]), alt=0, body=None):
        return GeocentricState(self, lon, lat, alt, velocity, body, epoch)
    
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

