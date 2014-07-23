import math
import numpy

from . import time

def deg_to_rad(deg):
    return (deg / 180.0)*math.pi

def rad_to_deg(rad):
    return (rad / math.pi)*180.0

class LatLonAlt(object):
    def __init__(self, body, lat_rad, lon_rad, alt_m=0):
        self.body = body
        self.lat = lat_rad
        self.lon = lon_rad
        self.alt = alt_m
    
    def to_xyz(self, date):
        day_f = time.seconds(date / self.body.day_length)
        lat = self.lat
        lon = self.lon + 2*math.pi*day_f
        alt = self.alt + self.body.eq_radius
        x_comp = math.sin(lon) * math.cos(lat)
        y_comp = math.cos(lon) * math.cos(lat)
        z_comp = math.sin(lat)
        return np.array([x_comp*alt, y_comp*alt, z_comp*alt])
