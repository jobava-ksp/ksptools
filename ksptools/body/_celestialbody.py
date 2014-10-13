import pickle
from numpy import zeros
from ._body import Body
from ._atmosphere import parse_atmosphere
from ._staticlocal import StaticSite
from .._kepler import KeplerOrbit as kepler
from .._math import asunits
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
        if hasattr(frame, 'orbit'):
            self.period = frame.orbit.period
            self.rap = frame.orbit.rap
            self.rpe = frame.orbit.rpe
        else:
            self.period = None
            self.rap = 0
            self.rpe = 0
    
    def statevector(self, t):
        return self.frame.orbit.statevector_by_time(t)
    
    def atmstate_by_statevector(self, stv, t):
        if self.atmosphere is None:
            return 0, 0, zeros(3), zeros(3)
        else:
            return self.atmosphere.atmstate_by_statevector(stv, t)
    
    def atmstate_by_lla(self, lat, lon, alt, t):
        if self.atmosphere is None:
            return 0, 0, zeros(3), zeros(3)
        else:
            return self.atmosphere.atmstate_by_lla(lat, lon, alt, t)
    
    def llav(self, stv, t):
        return self.surface.geodetic_llav(stv, t)
    
    def surface_uniti(self, lat, lon, t):
        return self.surface.uniti(lat, lon, t)
    
    def surface_unitj(self, lat, lon, t):
        return self.surface.unitj(lat, lon, t)
    
    def surface_unitk(self, lat, lon, t):
        return self.surface.unitk(lat, lon, t)
    
    def ijk_by_ll(self, lat, lon, t):
        return self.surface.uniti(lat,lon,t), self.surface.unitj(lat,lon,t), self.surface.unitk(lat,lon,t)
    
    def synodic(self, other):
        return kepler.synodic(self.frame.orbit, other.frame.orbit)


class System(object):
    def __init__(self):
        self.nodes = dict()
        self.keynames = dict()
        self.names = dict()
    
    def _addnode(self, keyname, name, n):
        self.nodes[keyname] = n
        self.keynames[name] = keyname
        self.names[keyname] = name
    
    def sun(self, keyname, name, u, geodetic_expr):
        cbody = CelestialBody(None, inertial_frame(), u, float('+inf'), parse_geodetic_frame(geodetic_expr))
        self._addnode(keyname, name, cbody)
        return cbody
    
    def cbody(self, keyname, name, parent_keyname, u, soi, geodetic_expr, orbit_expr, atm_expr=None):
        cbody = CelestialBody(
            self.nodes[parent_keyname],
            parse_orbital_frame(orbit_expr, self.nodes[parent_keyname].GM, 0),
            u,
            soi,
            parse_geodetic_frame(geodetic_expr))
        if atm_expr is not None:
            cbody.atmosphere = parse_atmosphere(cbody, atm_expr)
        self._addnode(keyname, name, cbody)
        return cbody
    
    def site(self, keyname, name, parent_keyname, lla_expr):
        lat, lon, alt = asunits(lla_expr[1:-1].split(','),['rad', 'rad', 'm'])
        site = StaticSite(self.nodes[parent_keyname], lat, lon, alt)
        self._addnode(keyname, name, site)
        return site
    
    def save(self, fname):
        with open(fname, 'w') as f:
            pickle.dump(self, f)
    
    def __getitem__(self, keyname):
        return self.nodes[keyname]
    
    @staticmethod
    def load(fname):
        with open(fname, 'r') as f:
            return pickle.load(f)
    

