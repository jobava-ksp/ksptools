from numpy import array, cos, sin, arctan
from numpy.linalg import norm

from . import orbit
from .util import Ax, rotvec, unitz, unity, unitz

class State(object):
    def __init__(self, refbody, body=None, epoch=0.):
        self.refbody = refbody
        self.body = body
        self.epoch = epoch
        self.up = unitz #TODO: more generalized body orientation
    
    def asorbit(self):
        raise NotImplementedError
    
    def asvectors(self):
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
                self.refbody.std_g_param,
                self.epoch)
        return OrbitalState(self.refbody, kep, self.body, self.epoch)

    def asvectors(self):
        return self
    
    def asgeocentric(self):
        # rotation
        up = self.refbody.state.up
        t = (2*pi/(self.refbody.sidereal_rate)) * self.epoch
        Ar = rotvec(u,-t)
        r = Ax(Ar, self.position)
        
        # fixed position
        alt = norm(r) - self.refbody.eq_radius
        x,y,z = r
        lat = arcsin(z/norm(r))
        lon = arccos(y/r)
        if y < 0:
            lon = 2*pi - lon
        
        # rotate velocity
        Av = rotvec(unity, -lat)
        v = Ax(Av*Ar, self.velocity)
        
        return GeocentricState(self.refbody, lon, lat, alt, v, self.body, self.epoch)
    
    def rv(self, t):
        return self.position + self.velocity * (t-epoch), self.velocity


class GeocentricState(State):
    def __init__(self, refbody, longitude, latitude, altitude, velocity, body=None, epoch=0.):
        State.__init__(self, refbody, body, epoch):
        self.longitude = longitude
        self.latitude = latitude
        self.altitude = altitude
        self.velocity = velocity
    
    def asvectors(self):
        r = self.altitude + self.refbody.eq_radius
        v = self.velocity
        
        # fixed position
        x = r*cos(self.latitude)*cos(self.longitude)
        y = r*cos(self.latitude)*sin(self.longitude)
        z = r*sin(self.latitude)
        fixed = array([x,y,z])
        
        # rotation
        up = self.refbody.state.up
        t = (2*pi/(self.refbody.sidereal_rate)) * self.epoch
        Ar = rotvec(u,t)
        Av = rotvec(unity,self.latitude)
        
        return VectorState(self.refbody, Ax(Ar, fixed), Ax(Ar*Av, v), self.body, self.epoch)
    
    def asgeocentric(self):
        return self


class FixedState(State):
    def __init__(self, refbody=None, body=None, epoch=0.):
        State.__init__(self, refbody, body, epoch)
