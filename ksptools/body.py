import re

from . import orbit
from . import atmosphere


class Body(object):
    def __init__(self, name, body_orbit, mass=None, u=None):
        self.name = name
        if u is None:
            self.std_g_param = mass * 6.67384e-11
            self.mass = mass
        else:
            self.std_g_param = u
            self.mass = u / 6.67384e-11
        self.orbit = body_orbit


class CelestialBody(Body):
    def __init__(self, name, eq_radius, u, sidereal_rate, soi, body_atmosphere, body_orbit):
        Body.__init__(self, name, body_orbit, u=u)
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
        if conf_parser.has_section(orbit_section_name):
            body_orbit = orbit.KeplerOrbit.from_config(conf_parser, orbit_section_name, parent)
        else:
            body_orbit = None
        if conf_parser.has_section(atm_section_name):
            body_atmosphere = atmosphere.Atmosphere.from_config(conf_parser, atm_section_name)
        else:
            body_atmosphere = None
        name = conf_parser.get(section, 'name')
        eq_radius = conf_parser.getfloat(section, 'radius')
        u = conf_parser.getfloat(section, 'gravitational_param')
        sidereal_rate = conf_parser.getfloat(section, 'sidereal_rate')
        soi = conf_parser.getfloat(section, 'soi')
        return key_name, cls(name, eq_radius, u, sidereal_rate, soi, body_atmosphere, body_orbit)
        
        
