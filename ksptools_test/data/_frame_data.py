import itertools

from numpy import arange, array, cos, dot, pi, sin, sqrt

from ksptools._frame import surface_frame, geodetic_frame
from ksptools._math import uniti, unitj, unitk, roty, rotz


#             #
# Surface ijk #
#             #

def obj_surface_ijk_func(s, t):
    i,j,k = s.uniti(t), s.unitj(t), s.unitk(t)
    return array([i, j, k])

testset_surface_ijk = list()

def make_surface_ijk_by_ll(lat, lon):
    g = geodetic_frame(3, 3, 0, 0, 0, 4.0)
    jki = array(map(lambda x: dot(roty(-lat), x).A1, (unitj, unitk, uniti)))
    s = surface_frame(g, lat, lon, 1)
    
    ffunc = lambda x, t: dot(rotz(t * (pi/2) + lon),x).A1
    gfunc = lambda t: lambda x: ffunc(x, t)
    hfunc = lambda s, t: ((s, t), array(map(gfunc(t), jki)))

    return list(hfunc(s, t) for t in arange(0, 4, 0.25))

for lat, lon in itertools.product((pi/4, 0, -pi/4), arange(0, 2*pi, pi/4)):
    testset_surface_ijk += make_surface_ijk_by_ll(lat, lon)

#                         #
# Surface inertial vector #
#                         #

def obj_surface_inertial_statevector_func(g, lat, lon, alt, t):
    return array(g.surface_inertial_statevector(lat, lon, alt, t).rv)

testset_surface_inertial_statevector = list()

def make_surface_inertial_statevector(g, lat, lon):
    t_alt = list(itertools.product(arange(0, 4, 0.25), [0,1]))
    
    def xfunc(alt, t):
        return (g, lat, lon, alt, t)
    
    def yfunc(alt, t):
        st = g._w*t + lon
        r = (g.surface_height(lat)+alt)*array([cos(st)*cos(lat), sin(st)*cos(lat), sin(lat)])
        v = (g._w*(g.surface_height(lat)+alt))*array([-sin(st)*cos(lat), cos(st)*cos(lat), 0])
        return array([r,v])
        
    return list((xfunc(alt, t), yfunc(alt, t)) for t, alt in t_alt)

g = geodetic_frame(3, 3, 0, 0, 0, 4.0)
for lat, lon in itertools.product((pi/4, 0, -pi/4),arange(0,2*pi,pi/4)):
    testset_surface_inertial_statevector += make_surface_inertial_statevector(g, lat, lon)

#               #
# Goedetic llav #
#               #

def obj_geodetic_llav_func(g, stv, t):
    lat, lon, alt, _ = g.geodetic_llav(stv, t)
    return array([lat, lon, alt])

def make_geodetic_llav(g):
    lat_range = [pi/4, 0, -pi/4]
    lon_range = list(arange(0, 2*pi, 0.25*pi))
    alt_range = [0, 1]
    t_range = list(arange(0, 4, 0.25))
    
    x_range = list(itertools.product(lat_range, lon_range, alt_range, t_range))
    
    x = [(g, g.surface_inertial_statevector(*t), t[-1]) for t in x_range]
    y = [array(t[0:3]) for t in x_range]
    return zip(x, y)

testset_geodetic_llav = make_geodetic_llav(geodetic_frame(3, 3, 0, 0, 0, 4.0))

