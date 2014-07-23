from . import orbit
from . import atmosphere

class Body(object):
    def __init__(self, name, mass):
        self.name = name
        self.mass = mass


class FixedBody(Body):
    def __init__(self, name, mass, body_orbit):
        Body.__init__(self, name, mass)
        self.orbit = orbit


class CelestialBody(FixedBody):
    def __init__(self, name, eq_radius, surface_area, mass, std_g_param, density, surface_g, escape_v, sidereal_rate, day_length, sidereal_v, soi, body_atmosphere, body_orbit):
        FixedBody.__init__(self, name, mass, body_orbit)
        self.eq_radius = eq_radius
        self.surface_area = surface_area
        self.std_g_param = std_g_param
        self.density = density
        self.surface_g = surface_g
        self.escape_v = escape_v
        self.sidereal_rate = sidereal_rate
        self.sidereal_v = sidereal_v
        self.day_length = day_length
        self.soi = soi
        self.atmosphere = body_atmosphere
