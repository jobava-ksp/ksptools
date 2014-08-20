from numpy import array, cos, sin, arccos, arcsin, arctan
from numpy.linalg import norm

from . import orbit
from .util import Ax, rotvec, uniti, unitj, unitk


def rv_to_llav(pos, vel, refbody, epoch):
        alt = norm(pos) - refbody.eq_radius
        t = ((2*pi/refbody.sidereal_rate) * epoch) % (2*pi)  # angle of planetary rotation
        R = rotz(t)                                          # rotation matrix for planetary rotation
        
        r = Ax(R.T, pos)                                     # rotate position into local position at t0
        v = Ax(R.T, vel)                                     # rotate velocity into local velocity at t0
        
        lon = arccos(r[0]/norm(r))
        lat = arcsin(r[2]/norm(r))
        if lat < 0:
            lon = -lon
        
        L = rotz(lon)*roty(-lat)                             # rotation matrix for gecentric -> local position (at t0)
        v = Ax(L.T, v)
        
        return lon, lat, alt, v

def llav_to_rv(lon, lat, alt, vel, refbody, epoch):
        rad = alt + refbody.eq_radius
        t = ((2*pi/refbody.sidereal_rate) * epoch) % (2*pi)  # angle of planetary rotation
        R = rotz(t)                                          # rotation matrix for planetary rotation
        L = rotz(lon)*roty(-lat)                             # rotation matrix for gecentric -> local position (at t0)
        
        r = Ax(R*L, rad*uniti)
        v = Ax(R*L, vel)
        return r, v

class State(object):
    def __init__(self, refbody, body=None, epoch=0.):
        self.refbody = refbody
        self.body = body
        self.epoch = epoch
    
    def _asorbit(self):
        raise NotImplementedError
    
    def _asvectors(self):
        raise NotImplementedError
    
    def _asgeocentric(self):
        raise NotImplementedError
    
    @classmethod
    def _from_vectors(cls, refbody, r, v, body, epoch):
        raise NotImplementedError
    
    @property
    def apoapsis_hieght(self):
        return norm(self.apoapsis)
    
    @property
    def apoapsis_altitude(self):
        return norm(self.apoapsis) - self.refbody.eq_radius
    
    @property
    def periapsis_height(self):
        return norm(self.periapsis)
    
    @property
    def periapsis_altitude(self):
        return norm(self.periapsis) - self.refbody.eq_radius
    
    @property
    def time_to_apoapsis(self):
        return self._asorbit().time_to_apoapsis
    
    @property
    def time_to_periapsis(self):
        return self._asorbit().time_to_periapsis
    
    @property
    def apoapsis(self):
        return self._asorbit().apoapsis
    
    @property
    def periapsis(self):
        return self._asorbit().periapsis
    
    @property
    def longitude(self):
        return self._asgeocentric().longitude
    
    @property
    def latitude(self):
        return self._asgeocentric().latitude
    
    @property
    def altitude(self):
        return norm(self.position) - self.refbody.eq_radius
    
    @property
    def kepler(self):
        return self._asorbit().kepler
    
    @property
    def period(self):
        return self._asorbit().kepler.period()


class OrbitalState(State):
    def __init__(self, refbody, kepler, body=None, epoch=0.):
        State.__init__(self, refbody, body, epoch)
        self._kepler = kepler
    
    def _asorbit(self):
        return self
    
    def _asvectors(self):
        if not hasattr(self, '_vectors'):
            r, v = self._kepler.rv(self.epoch)
            self._vectors = VectorState(self.refbody, r, v, self.body, self.epoch)
        return self._vectors
    
    def _asgeocentric(self):
        if not hasattr(self, '_geocentric'):
            r, v = self._kepler.rv(self.epoch)
            lon, lat, alt, vel = rv_to_llav(r, v, self.refbody, self.epoch)
            self._geocentric = GeocentricState(self.refbody, lon, lat, alt, vel, self.body, self.epoch)
        return self._geocentric
    
    @classmethod
    def from_rv(cls, refbody, r, v, body, epoch):
        kep = KeplerOrbit.from_rvu(r, v, refbody.std_g_param, epoch)
        return cls(refbody, kep, body, epoch)
    
    @property
    def position(self):
        return self._kepler.r(self.epoch)
    
    @property
    def velocity(self):
        return self._kepler.v(self.epoch)
    
    @property
    def rv(self):
        return self._kepler.rv(self.epoch)
    
    @property
    def apoapsis(self):
        return self._kepler.ap()
    
    @property
    def periapsis(self):
        return self._kepler.pe()
    
    @property
    def time_to_apoapsis(self):
        return self._kepler.time_to_ap(self.epoch)
    
    @property
    def time_to_periapsis(self):
        return self._kepler.time_to_pe(self.epoch)
    
    @property
    def kepler(self):
        return self._kepler
    
    @property
    def period(self):
        return self._kepler.period()


class VectorState(State):
    def __init__(self, refbody, r, v, body=None, epoch=0.):
        State.__init__(self, refbody, body, epoch)
        self._vel = v
        self._pos = r
    
    def _asorbit(self):
        if not hasattr(self, '_orbit'):
            kep = orbit.KeplerOrbit.from_rvu(
                    self._pos,
                    self._vel,
                    self.refbody.std_g_param,
                    self.epoch)
            self._orbit = OrbitalState(self.refbody, kep, self.body, self.epoch)
        return self._orbit

    def _asvectors(self):
        return self
    
    def _asgeocentric(self):
        if not hasattr(self, '_gocentric'):
            lon, lat, alt, vel = rv_to_llav(self._pos, self._vel, self.refbody, self.epoch)
            self._geocentric = GeocentricState(self.refbody, lon, lat, alt, vel, self.body, self.epoch)
        return self._geocentric
    
    @classmethod
    def from_rv(cls, refbody, r, v, body, epoch):
        return cls(refbody, r, v, body, epoch)
    
    @property
    def position(self):
        return self._pos
    
    @property
    def velocity(self):
        return self._vel
    
    @property
    def rv(self):
        return self._pos, self._vel


class GeocentricState(State):
    def __init__(self, refbody, longitude, latitude, altitude, velocity, body=None, epoch=0.):
        State.__init__(self, refbody, body, epoch)
        self._lon = longitude
        self._lat = latitude
        self._alt = altitude
        self._vel = velocity
    
    def _asorbit(self):
        if not hasattr(self, '_orbit'):
            r, v = llav_to_rv(self._lon, self._lat, self._alt, self._vel, self.refbody, self.epoch)
            kep = orbit.KeplerOrbit.from_rvu(r, v, self.refbody.std_g_param, self.epoch)
            self._orbit = OrbitalState(self.refbody, kep, self.body, self.epoch)
    
    def _asvectors(self):
        r, v = llav_to_rv(self._lon, self._lat, self._alt, self._vel, self.refbody, self.epoch)
        return VectorState(self.refbody, Ax(Ar, fixed), Ax(Ar*Av, v), self.body, self.epoch)
    
    def _asgeocentric(self):
        return self
    
    @classmethod
    def from_rv(cls, refbody, r, v, body, epoch):
        lon, lat, alt, vel = rv_to_llav(r, v, refbody, refbody, refbody.pu, epoch)
        return cls(self.refbody, lon, lat, alt, vel, body, epoch)
    
    @property
    def position(self):
        r, v = llav_to_rv(self._lon, self._lat, self._alt, self._vel, self.refbody, self.epoch)
        return r
    
    @property
    def velocity(self):
        r, v = llav_to_rv(self._lon, self._lat, self._alt, self._vel, self.refbody, self.epoch)
        return v
    
    @property
    def rv(self):
        r, v = llav_to_rv(self._lon, self._lat, self._alt, self._vel, self.refbody, self.epoch)
        return r, v


class FixedState(State):
    def __init__(self, refbody=None, body=None, epoch=0.):
        State.__init__(self, refbody, body, epoch)
    
    @property
    def position(self):
        return array([0,0,0])
    
    @property
    def velocity(self):
        return array([0,0,0])
