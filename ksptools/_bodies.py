from . import body
from . import orbit

def _orbit(parent, a, e, i, argpe, la, M, epoch):
    return orbit.KeplerOrbit.from_planet_paremters(parent, a, e, i, argpe, la, M, epoch)

def _atm(psl, atmsl, scale_height, height, temp_min, temp_max, hasO2):
    return atmposphere.Atmposhere(psl, atmsl, scale_height, height, temp_min, temp_max, hasO2)

def _planet(name, parent, eq_radius, surface_A, mass, std_g_param, density, surface_g, escape_v, sidereal_rate, day_length, sidereal_v, soi, _p_orbit=None, _p_atm=None):
    return body.CelestialBody(name, eq_radius, surface_A, mass, std_g_param, density, surface_g, escape_v, sidereal_rate, day_length, sidereal_v, soi, _p_atm, _p_orbit)

_kerbol = _planet(
    'Kerbol',
    None,
    2.616*10**8,
    8.5997404*10**17,
    1.7565670*10**28,
    1.1723328*10**18,
    234.24098,
    17.1,
    94672.01,
    4.32*10**5,
    0,
    3804.8,
    float('inf'))

_moho = _planet(
    'Moho',
    
