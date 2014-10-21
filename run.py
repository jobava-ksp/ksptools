from ksptools import *


import numpy as np
import numpy.linalg as la


poodle = partialstage('poodle', 2.0e+3, 220e+3, 390, 270, tankstep=4.0e+3)
skipper = partialstage('skipper', 3.0e+3, 650e+3, 370, 320, tankstep=4.0e+3)
mainsail = partialstage('mainsail', 6.0e+3, 1500e+3, 360, 320, tankstep=4.0e+3)

bacc = partialbooster('BACC', 1.505e+3, 315e+3, 250, 230, tankstep=0, tmratio=0)
kd25k = partialbooster('S1 SRB-KD25k', 3.0e+3, 650e+3, 250, 230, tankstep=0, tmratio=0)


for booster_type, booster_fuel in [(bacc, 6.37e+3), (kd25k, 18.75e+3)]:
    for n in [0,2,3,4,6]:
        if n > 0:
            print('\t{} {} boosters:'.format(n, booster_type.name))
            twr = [1.5, 1.7, 2.0]
            stages = [skipper, mainsail, (n * booster_type)]
            booster_stage_fuel = n * booster_fuel
            second_stage_fuel = booster_type.burntime(booster_fuel)*mainsail.ff()
            print('\tStage Fuel: {:.3e}kg'.format(booster_stage_fuel))
            fixed = [(1, second_stage_fuel, None),
                     (2, booster_stage_fuel, booster_stage_fuel)]
        else:
            print('\tNo boosters:')
            twr = [0.5, 1.7, 1.7]
            stages = [poodle, skipper, mainsail]
            fixed = []
        for name, mp, stwr, dv, sumdv in minimizefuel(10.5e+3, stages, twr, 6200, fixed=fixed):
            print('\t\t{}: {:.6e}kg fuel, {:.2f} TWR, {:.1f}m/s | {:.1f}m/s'.format(name,mp,stwr,dv,sumdv))

