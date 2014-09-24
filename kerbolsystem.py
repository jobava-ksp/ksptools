
#
ks = System()
ks.sun('kerbol', 'Kerbol', '<2.616e+8m, 2.616e+8m, 0, 0, 0, 4.32e+5sec>', 1.1723328e+18, float('+inf'))

ks.cbody('moho', 'Moho', 'kerbol', '<2.5e+5m, 2.5e+5m, 0, 0, 0, 1.21e+6sec>', 1.6860938e+11, 9.646663e+6)
ks.cbody('eve', 'Eve', 'kerbol', '<7.0e+5m, 7.0e+5m, 0, 0, 0, 8.05e+5sec>', 8.1717302e+12, 8.5109365e+7)
ks.cbody('gilly', 'Gilly', 'eve', '<1.3e+4m, 1.3e+4m, 0, 0, 0, 2.8255e+4sec>', 8.2894498e+6, 1.2612327e+5)
ks.cbody('kerbin', 'Kerbin', 'kerbol', '<6.0e+5m, 6.0e+5m, 0, 0, 0, 2.16e+4sec>', 3.5316e+12, 8.4169286e+7)
ks.cbody('mun', 'Mun', 'kerbin', '<2.0e+5m, 2.0e+5m, 0, 0, 0, 1.38984e+5sec>', 6.5138398e+10, 2.4295591e+6)
ks.cbody('minmus', 'Minmus', 'kerbin', 
#
