import math

class Atmposhere(object):
    def __init__(self, pressure_at_sl, atm_at_sl, scale_height, height, temp_min, temp_max, hasO2):
        self.pressure_at_sl = pressure_at_sl
        self.scale_height = scale_height
        self.height = height
        self.temp_min = temp_min
        self.temp_max = temp_max
        self.hasO2 = hasO2
    
    def atm(self, altitude):
        return self.atm_at_sl * (math.e ** (-altitude / self.scale_height))
    
    def p(self, altitude):
        return self.pressure_at_sl * (math.e ** (-altitude / self.scale_height))

