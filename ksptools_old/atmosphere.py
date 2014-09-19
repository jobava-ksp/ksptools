import math

from . import unit

class Atmosphere(object):
    def __init__(self, pressure_at_sl, atm_at_sl, scale_height, height, temp_min, temp_max, hasO2):
        self.pressure_at_sl = pressure_at_sl
        self.atm_at_sl = atm_at_sl
        self.scale_height = scale_height
        self.height = height
        self.temp_min = temp_min
        self.temp_max = temp_max
        self.hasO2 = hasO2
    
    def atm(self, altitude):
        return self.atm_at_sl * (math.e ** (-altitude / self.scale_height))
    
    def p(self, altitude):
        return self.pressure_at_sl * (math.e ** (-altitude / self.scale_height))
    
    @classmethod
    def from_config(cls, config_parser, section):
        psl = config_parser.getfloat(section, 'psl')
        asl = config_parser.getfloat(section, 'asl')
        scale_height = config_parser.getfloat(section, 'scale_height')
        height = config_parser.getfloat(section, 'height')
        temp_min = unit.tounit(config_parser.get(section, 'temp_min'), 'C')
        temp_max = unit.tounit(config_parser.get(section, 'temp_max'), 'C')
        hasO2 = config_parser.getboolean(section, 'has_oxygen')
        return cls(psl,asl,scale_height,height,temp_min,temp_max,hasO2)

