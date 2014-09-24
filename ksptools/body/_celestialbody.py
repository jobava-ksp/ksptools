from ._body import Body
from .._math import parse_kepler
from .._frame import geocentric_frame, geodetic_frame, orbital_frame, inertial_frame

class CelestialBody(Body):
    def __init__(self, parent_node, frame, u, eq_radius, polar_radius, inc, lonasc, argve, sidereal_period):
        Body.__init__(parent_node, frame)
        self.u = u
        self.equitorial_radius = eq_radius
        self.geocentric_frame = geocentric_frame(inc, lonasc, argve)
        self.surface_frame = geodetic_frame(polar_radius, eq_radius, inc, lonasc, argve, sidereal_period)

