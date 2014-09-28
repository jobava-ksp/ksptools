import numpy as np
from .._math import asunits
from ._field import Field


class Atmosphere(Field):
    def __init__(self, parent_node, p_sl, atm_sl, height, scale_height, min_temp, max_temp, has_oxygen):
        Field.__init__(self, parent_node)
        self.p_sl = p_sl
        self.atm_sl = atm_sl
        self.height = height
        self.scale_height = scale_height
        self.min_temp = min_temp
        self.max_temp = max_temp
        self.has_oxygen = has_oxygen
        self.surface = parent_node.surface
    
    def atm_state(self, stv, t):
        lat, lon, alt, v = self.surface.geodetic_llav(stv, t)
        if alt < self.height:
            p = self.p_sl * pow(np.e, -alt/self.scale_height)
            a = self.atm_sl * pow(np.e, -alt/self.scale_height)
            va = self.surface.surface_inertial_statevector(lat, lon, alt, t)
            return p, a, v, stv.v - va
        else:
            return 0, 0, np.zeros(3), np.zeros(3)
    

def parse_atmosphere(parent, atm_expr):
    params = [e.strip() for e in atm_expr[1:-1].split(',')]
    params = [parent] + list(map(float,params[0:2])) + list(asunits(params[2:6], ['m','m','C','C'])) + [bool(params[6])]
    return Atmosphere(*params)
    
