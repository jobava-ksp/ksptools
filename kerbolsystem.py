from ksptools import system

def make_system():
    ks = system()
    kerbol = ks.sun(
            'kerbol', 'Kerbol', 1.1723328e+18,
            '<2.616e+8m, 2.616e+8m, 0, 0, 0, 4.32e+5sec>')

    moho = ks.cbody(
            'moho', 'Moho', 'kerbol', 1.6860938e+11, 9.646663e+6,
            '<2.5e+5m, 2.5e+5m, 0, 0, 0, 1.21e+6sec>',
            '<5.263138304e+9m, 0.2, 7.0deg, 15.0deg, 70.0deg, 3.14rad>')
            
    eve = ks.cbody(
            'eve', 'Eve', 'kerbol', 8.1717302e+12, 8.5109365e+7,
            '<7.0e+5m, 7.0e+5m, 0, 0, 0, 8.05e+5sec>',
            '<9.832684544e+9m, 0.01, 2.1deg, 0.0deg, 15.0deg, 3.14rad>',
            '<5.06625e+2, 5.0, 7.0e+3m, 9.6708574e+4m, -40.19C, 149.96C, False>')
            
    gilly = ks.cbody(
            'gilly', 'Gilly', 'eve', 8.2894498e+6, 1.2612327e+5,
            '<1.3e+4m, 1.3e+4m, 0, 0, 0, 2.8255e+4sec>',
            '<3.15e+7m, 0.55, 12.0deg, 10.0deg, 80deg, 0.9rad>')
            
    kerbin = ks.cbody(
            'kerbin', 'Kerbin', 'kerbol', 3.5316e+12, 8.4169286e+7,
            '<6.0e+5m, 6.0e+5m, 0, 0, 0, 2.16e+4sec>',
            '<1.3599840256e+10m, 0, 0, 0, 0, 3.14rad>',
            '<1.01325e+2, 1.0, 5.0e+3m, 6.9077553e+4m, -40.19C, 20.0C, True>')
            
    mun = ks.cbody(
            'mun', 'Mun', 'kerbin', 6.5138398e+10, 2.4295591e+6,
            '<2.0e+5m, 2.0e+5m, 0, 0, 0, 1.38984e+5sec>',
            '<1.2e+7m, 0, 0, 0, 0, 1.7rad>')
            
    minmus = ks.cbody(
            'minmus', 'Minmus', 'kerbin', 1.7658e+9, 2.2474284e+5,
            '<6.0e+4m, 6.0e+4m, 0, 0, 0, 4.04e+4sec>',
            '<4.7e+7m, 0, 6.0deg, 38.0deg, 78.0deg, 0.9rad>')
            
    duna = ks.cbody(
            'duna', 'Duna', 'kerbol', 3.0136321e+11, 4.7921949e+7,
            '<3.2e+5m, 3.2e+5m, 0, 0, 0, 6.5517859e+4sec>',
            '<2.0726155264e+10m, 0.05, 0.06deg, 0, 135.5deg, 3.14rad>',
            '<2.02650e+1, 0.2, 3.0e+3m, 4.1446532e+4m, -50.24C, -30.17C, False>')
            
    ike = ks.cbody(
            'ike', 'Ike', 'duna', 1.8568369e+10, 1.0495989e+6,
            '<1.3e+5m, 1.3e+5m, 0, 0, 0, 6.5517862e+4sec>',
            '<3.2e+6m, 0.03, 0.2deg, 0.0deg, 0.0deg, 1.7rad>')
            
    dres = ks.cbody(
            'dres', 'Dres', 'kerbol', 2.1484489e+10, 3.2832840e+7,
            '<1.38e+5m, 1.38e+5m, 0, 0, 0, 3.48e+4sec>',
            '<4.0839348203e+10m, 0.14, 5.0deg, 90deg, 280deg, 3.14rad>')
            
    jool = ks.cbody(
            'jool', 'Jool', 'kerbol', 2.82528e+14, 2.4559852e+9,
            '<6.0e+6m, 6.0e+6m, 0, 0, 0, 3.6e+4sec>',
            '<6.8773560320e+10m, 0.05, 1.304deg, 0deg, 52deg, 0.1rad>',
            '<1.51988e+4, 15.0, 1.0e+4m, 1.3815511e+5m, -86.13C, 976.55C, False>')
            
    laythe = ks.cbody(
            'laythe', 'Laythe', 'jool', 1.962e+12, 3.7236458e+6,
            '<5.0e+5m, 5.0e+5m, 0, 0, 0, 5.2980879e+4sec>',
            '<2.7184e+7m, 0, 0, 0, 0, 3.14rad>',
            '<8.106e+1, 0.8, 4.0e+3m, 5.5262042e+4m, -40.19C, 6.21C, True>')
            
    vall = ks.cbody(
            'vall', 'Vall', 'jool', 2.074815e+11, 2.4064014e+6,
            '<3.0e+5m, 3.0e+5m, 0, 0, 0, 1.0596209e+5sec>',
            '<4.3152e+7m, 0, 0, 0, 0, 0.9rad>')
            
    tylo = ks.cbody(
            'tylo', 'Tylo', 'jool', 2.82528e+12, 1.0856518e+7,
            '<6.0e+5m, 6.0e+5m, 0, 0, 0, 2.1192636e+5sec>',
            '<6.85e+7m, 0, 0.025deg, 0, 0, 3.14deg>')
            
    bop = ks.cbody(
            'bop', 'Bop', 'jool', 2.4868349e+9, 1.2210609e+6,
            '<6.5e+4m, 6.5e+4m, 0, 0, 0, 5.445074e+5sec>',
            '<1.285e+8m, 0.24, 15.0deg, 25.0deg, 10.0deg, 0.9rad>')
            
    pol = ks.cbody(
            'pol', 'Pol', 'jool', 7.2170208e+8, 1.0421389e+6,
            '<4.4e+4m, 4.4e+4m, 0, 0, 0, 9.0190262e+5sec>',
            '<1.7989e+7m, 0.17, 4.25deg, 15.0deg, 2.0deg, 0.9rad>')

    eeloo = ks.cbody(
            'eeloo', 'Eeloo', 'kerbol', 7.4410815e+10, 1.1908294e+8,
            '<2.1e+5m, 2.1e+5m, 0, 0, 0, 1.9460e+4sec>',
            '<9.011882e+10m, 0.26, 6.15deg, 260deg, 50deg, 3.14rad>')

    #ks.export('kerbol_system.pickle')
    return ks

kerbolsystem = make_system()

