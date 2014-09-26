import pickle
from numpy import zeros
from ._body import Body
from ._atmosphere import parse_atmosphere
from .._vector import statevector
from .._frame import inertial_frame
from .._frame import parse_geodetic_frame, parse_orbital_frame


class CelestialBody(Body):
    def __init__(self, parent_node, frame, u, soi, geodetic_surface):
        Body.__init__(self, parent_node, frame)
        self.GM = u
        self.surface = geodetic_surface
        self.surface_frame = geodetic_surface.frame
        self.atmosphere = None
        self.soi = soi
    
    def statevector(self, t):
        return self.frame.toinertial(statevector(zeros(3), zeros(3)), t)


class System(object):
    def __init__(self):
        self.bodies = dict()
        self.keynames = dict()
        self.names = dict()
    
    def _addbody(self, keyname, name, body):
        self.bodies[keyname] = body
        self.keynames[name] = keyname
        self.names[keyname] = name
    
    def sun(self, keyname, name, u, geodetic_expr):
        cbody = CelestialBody(None, inertial_frame(), u, float('+inf'), parse_geodetic_frame(geodetic_expr))
        self._addbody(keyname, name, cbody)
        return cbody
    
    def cbody(self, keyname, name, parent_keyname, u, soi, geodetic_expr, orbit_expr, atm_expr=None):
        cbody = CelestialBody(
            self.bodies[parent_keyname],
            parse_orbital_frame(orbit_expr, u, 0),
            u,
            soi,
            parse_geodetic_frame(geodetic_expr))
        if atm_expr is not None:
            cbody.atmosphere = parse_atmosphere(cbody, atm_expr)
        self._addbody(keyname, name, cbody)
        return cbody
    
    def export(self, fname):
        with open(fname, 'w') as f:
            pickle.dump(self, f)
    
    def __getitem__(self, keyname):
        return self.bodies[keyname]
    
    @staticmethod
    def load(fname):
        with open(fname, 'r') as f:
            return pickle.load(f)
    

