from ksptools.part import Catalogue, ResourceType, CrewPodPartType, ProbePodPartType, FuelTankPartType, EnginePartType
from ksptools.part import FeedPartType

kerbalparts = Catalogue()

## ResourceType(name, title, unitcost, unitdensity, flowmode) ##

kerbalparts.add(ResourceType('charge', 'Electric Charge', 0, 0, 'universal'))
kerbalparts.add(ResourceType('liquid_fuel', 'Liquid Fuel', 0.8, 5, 'stage'))
kerbalparts.add(ResourceType('oxidizer', 'Oxidizer', 0.18, 5, 'stage'))
kerbalparts.add(ResourceType('intake_air', 'Intake Air', 0, 5, 'universal'))
kerbalparts.add(ResourceType('solid_fuel', 'Solid Fuel', 0.6, 7.5, 'part'))
kerbalparts.add(ResourceType('monoprop', 'Mono Propellant', 1.2, 4, 'priority'))
kerbalparts.add(ResourceType('eva', 'EVA Propellant', 0, 0, 'part'))
kerbalparts.add(ResourceType('xenon', 'Xenon Gas', 4, 0.1, 'priority'))

## CrewPodPartType(name, title, radialsize, drycost, drymass, ceofdrag, kerbals, torque, resources) ##

kerbalparts.add(CrewPodPartType('mk1_cockpit', 'Mk1 Cockpit', [None, 'small'], 2191, 1.25e+3, 0.1, 1, 10, [kerbalparts.monoprop(7.5), kerbalparts.charge(50)]))
kerbalparts.add(CrewPodPartType('mk2_cockpit', 'Mk2 Cockpit', ['small', 'small'], 1591, 1.0e+3, 0.08, 1, 10, [kerbalparts.monoprop(7.5), kerbalparts.charge(50)]))
kerbalparts.add(CrewPodPartType('mk3_cockpit', 'Mk3 Cockpit', ['small', 'Mk3'], 3564, 3.5e+3, 0.15, 3, 10, [kerbalparts.monoprop(30.0), kerbalparts.charge(50)]))
kerbalparts.add(CrewPodPartType('mk1_pod', 'Command Pod Mk1', ['tiny', 'small'], 588.0, 0.8e+3, 0.15, 1, 5, [kerbalparts.monoprop(10.0), kerbalparts.charge(50)]))
kerbalparts.add(CrewPodPartType('mk12_pod', 'Mk1-2 Command Pod', ['small', 'large'], 3764, 4.0e+3, 0.15, 3, 15, [kerbalparts.monoprop(30.0), kerbalparts.charge(50)]))
kerbalparts.add(CrewPodPartType('mk1_lander', 'Mk1 Lander-can', ['small', 'small'], 1482, 0.6e+3, 0.2, 1, 3, [kerbalparts.monoprop(15.0), kerbalparts.charge(50)]))
kerbalparts.add(CrewPodPartType('mk2_lander', 'Mk2 Lander-can', ['large', 'large'], 3202, 2.5e+3, 0.15, 2, 15, [kerbalparts.monoprop(40.0), kerbalparts.charge(50)]))
kerbalparts.add(CrewPodPartType('copula', 'PPD-12 Copula Module', ['small', 'large'], 3188, 4.5e+3, 0.4, 1, 9, [kerbalparts.monoprop(10.0), kerbalparts.charge(200)]))
kerbalparts.add(CrewPodPartType('seat', 'EAS-1 External Command Seat', [None], 200, 0.05e+10, 0.05, 1, 0, []))

## ProbePodPartType(name, title, radialsize, drycost, drymass, ceofdrag, torque, resources) ##

kerbalparts.add(ProbePartType('qbe_probe', 'Probodobodyne QBE', ['tiny', 'tiny'], 600, 0.08e+3, 0.15, 0.5, [kerbalparts.charge(10)]))
kerbalparts.add(ProbePartType('hecs_probe', 'Probodobodyne HECS', ['tiny', 'tiny'], 450, 0.1e+3, 0.15, 0.5, [kerbalparts.charge(10)]))
kerbalparts.add(ProbePartType('okto_probe', 'Probodobodyne OKTO', ['tiny', 'tiny'], 450, 0.1e+3, 0.15, 0.3, [kerbalparts.charge(10)]))
kerbalparts.add(ProbePartType('okto2_probe', 'Probodobodyne OKTO2', ['tiny', 'tiny'], 480, 0.04e+3, 0.15, 0.5, [kerbalparts.charge(5)]))
kerbalparts.add(ProbePartType('stayput_probe', 'Stayputnik Mk. 1', [None, 'tiny'], 300, 0.05e+3, 0.15, 0.3, [kerbalparts.charge(10)]))
kerbalparts.add(ProbePartType('rc_small_probe', 'RC-001S Remote Guidance Unit', ['small', 'small'], 750, 0.1e+3, 0.2, 0.5, [kerbalparts.charge(30)]))
kerbalparts.add(ProbePartType('rc_large_probe', 'RC-L01 Remote Guidance Unit', ['large', 'large'], 2000, 0.5e+3, 0.2, 1.5, [kerbalparts.charge(30)]))

## FuelTankPartType(name, title, radialsize, drycost, drymass, coefdrag, resources) ##

kerbalparts.add(FuelTankPartType('mini_tank', 'Oscar-B Fuel Tank', ['tiny', 'tiny'], 174.15, 0.02e+3, 0.15, [kerbalparts.liquid_fuel(5.735), kerbalparts.oxidizer(7.0)]))
kerbalparts.add(FuelTankPartType('round_tank', 'ROUND-8 Toroidal Fuel Tank', ['small', 'small'], 349.8, 0.03e+10, 0.15, [kerbalparts.liquid_fuel(10.0), kerbalparts.oxidizer(12.2)]))
kerbalparts.add(FuelTankPartType('flt100_tank', 'FL-T100 Fuel Tank', ['small', 'small'], 204.1, 0.06e+3, 0.2, [kerbalparts.liquid_fuel(45), kerbalparts.oxidizer(55)]))
kerbalparts.add(FuelTankPartType('flt200_tank', 'FL-T200 Fuel Tank', ['small', 'small'], 333.2, 0.13e+3, 0.2, [kerbalparts.liquid_fuel(90), kerbalparts.oxidizer(110)]))
kerbalparts.add(FuelTankPartType('flt400_tank', 'FL-T400 Fuel Tank', ['small', 'small'], 666.4, 0.25e+3, 0.2, [kerbalparts.liquid_fuel(180), kerbalparts.oxidizer(220)]))
kerbalparts.add(FuelTankPartType('flt800_tank', 'FL-T800 Fuel Tank', ['small', 'small'], 1232.8, 0.5e+3, 0.2, [kerbalparts.liquid_fuel(360), kerbalparts.oxidizer(440)]))
kerbalparts.add(FuelTankPartType('x8_tank', 'Rockomax X200-8 Fuel Tank', ['large', 'large'], 1232.8, 0.5e+3, 0.2, [kerbalparts.liquid_fuel(360), kerbalparts.oxidizer(440)]))
kerbalparts.add(FuelTankPartType('x16_tank', 'Rockomax X200-16 Fuel Tank', ['large', 'large'], 2465.6, 1.0e+3, 0.2, [kerbalparts.liquid_fuel(720), kerbalparts.oxidizer(880)]))
kerbalparts.add(FuelTankPartType('x32_tank', 'Rockomax X200-32 Fuel Tank', ['large', 'large'], 4931.2, 2.0e+3, 0.2, [kerbalparts.liquid_fuel(1440), kerbalparts.oxidizer(1760)]))
kerbalparts.add(FuelTankPartType('x64_tank', 'Rockomax Jumbo-64 Fuel Tank', ['large', 'large'], 9862.4, 4.0e+3, 0.2, [kerbalparts.liquid_fuel(2880), kerbalparts.oxidizer(3520)]))
kerbalparts.add(FuelTankPartType('s36_tank', 'Kerbodyne S3-3600 Tank', ['exlarge', 'exlarge'], 5547.6, 2.5e+3, 0.2, [kerbalparts.liquid_fuel(1620), kerbalparts.oxidizer(1980)]))
kerbalparts.add(FuelTankPartType('s72_tank', 'Kerbodyne S3-7200 Tank', ['exlarge', 'exlarge'], 11095.2, 5.0e+3, 0.2, [kerbalparts.liquid_fuel(3240), kerbalparts.oxidizer(3960)]))
kerbalparts.add(FuelTankPartType('s144_tank', 'Kerbodyne S3-14400 Tank', ['exlarge', 'exlarge'], 16190.4, 10.0e+3, 0.2, [kerbalparts.liquid_fuel(6480), kerbalparts.oxidizer(7920)]))
kerbalparts.add(FuelTankPartType('strat5_cylinder_tank', 'Stratus-V Cylindrified Monopropellant Tank', [None], 620, 0.15e+3, 0.2, [kerbalparts.monoprop(150)]))
kerbalparts.add(FuelTankPartType('strat5_round_tank', 'Stratus-V Roundified Monopropellant Tank', [None], 352, 0.08e+3, 0.2, [kerbalparts.monoprop(40)]))
kerbalparts.add(FuelTankPartType('flr1', 'FL-R1 RCS Fuel Tank', ['large', 'large'], 400, 0.4e+3, 0.2, [kerbalparts.monoprop(750)]))
kerbalparts.add(FuelTankPartType('flr10', 'FL-R10 RCS Fuel Tank', ['tiny', 'tiny'], 340, 0.05e+3, 0.2, [kerbalparts.monoprop(50)]))
kerbalparts.add(FuelTankPartType('flr25', 'FL-R25 RCS Fuel Tank', ['small', 'small'], 680, 0.15e+3, 0.2, [kerbalparts.monoprop(100)]))

kerbalparts.add(FeedPartType('ftx', 'FTX-2 External Fuel Duct', [None], 150, 0.05e+3, 0.2))

## EnginePartType(name, title, radialsize, drycost, drymass, coefdrag, fueltypes, minthrust, maxthrust, isp, ispatm, resources=[]) ##

liquid_fuel_req = [(0.45, kerbalparts.liquid_fuel), (0.55, kerbalparts.oxidizer)]
solid_fuel_req = [(1.0, kerbalparts.solid_fuel)]

kerbalparts.add(EnginePartType('lv1', 'LV-1 Liquid Fuel Engine', ['tiny', 'tiny'], 350, 0.03e+3, 0.2, liquid_fuel_req, 0, 4e+3, 290, 220))
kerbalparts.add(EnginePartType('lv1r', 'LV-1R Liquid Fuel Engine', [None], 650, 0.03e+3, 0.2, liquid_fuel_req, 0, 4e+3, 290, 220))
kerbalparts.add(EnginePartType('rm24_77', 'Rockomax 24-77', [None], 480, 0.09e+3, 0.2, liquid_fuel_req, 0, 20e+3, 300, 250))
kerbalparts.add(EnginePartType('rm48_7s', 'Rockomax 48-7S', ['tiny', 'tiny'], 300, 0.1e+3, 0.2, liquid_fuel_req, 0, 30e+3, 350, 300))
kerbalparts.add(EnginePartType('rm55', 'Rockomax Mark 55 Radial Mount Liquid Engine', [None], 850, 0.9e+3, 0.2, liquid_fuel_req, 0, 120e+3, 320, 290))
## Aerospike engine ? ##
kerbalparts.add(EnginePartType('lvn', 'LV-N Atomic Rocket Motor', ['small', 'small'], 8700, 2.25e+3, 0.2, liquid_fuel_req, 0, 60e+3, 800, 220))
kerbalparts.add(EnginePartType('lv909', 'LV-909 Liquid Fuel Engine', ['small', 'small'], 750, 0.5e+3, 0.2, liquid_fuel_req, 0, 50e+3, 390, 300))
kerbalparts.add(EnginePartType('lv45', 'LV-45 Liquid Fuel Engine', ['small', 'small'], 950, 1.5e+3, 0.2, liquid_fuel_req, 0, 200e+3, 370, 320))
kerbalparts.add(EnginePartType('lv30', 'LV-30 Liquid Fuel Engine', ['small', 'small'], 850, 1.25e+3, 0.2, liquid_fuel_req, 0, 215e+3, 370, 320))
kerbalparts.add(EnginePartType('poodle', 'Rockomax \"Poodle\" Liquid Engine', ['large', 'large'], 1600, 2e+3, 0.2, liquid_fuel_req, 0, 220e+3, 390, 270))
kerbalparts.add(EnginePartType('skipper', 'Rockomax \"Skipper\" Liquid Engine', ['large', 'large'], 2850, 3e+3, 0.2, liquid_fuel_req, 0, 650e+3, 370, 320))
kerbalparts.add(EnginePartType('mainsail', 'Rockomax \"Mainsail\" Liquid Engine', ['large', 'large'], 5650, 6e+3, 0.2, liquid_fuel_req, 0, 1500e+3, 360, 320))
kerbalparts.add(EnginePartType('kr2l', 'Kerbodyne KR-2L Advanced Engine', ['exlarge', 'exlarge'], 20850, 6.5e+3, 0.2, liquid_fuel_req, 0, 2500e+3, 380, 280))
kerbalparts.add(EnginePartType('lfbx2', 'LFB KR-1x2', ['large', None], 13462, 10e+3, 0.2, liquid_fuel_req, 0, 2000e+3, 340, 290, [kerbalparts.liquid_fuel(2880), kerbalparts.oxidizer(3520)]))
kerbalparts.add(EnginePartType('ksx4', 'S3 KS-25x4 Engine Cluster', ['exlarge', None], 32400, 9.75e+3, 0.2, liquid_fuel_req, 0, 3200e+3, 360, 320))
kerbalparts.add(EnginePartType('sepratron', 'Sepratron I', [None], 45.2, 0.01e+3, 0.2, solid_fuel_req, 18e+3, 18e+3, 100, 100, [kerbalparts.solid_fuel(8)]))
kerbalparts.add(EnginePartType('rt10', 'RT-10 Solid Fuel Booster', ['small', 'small'], 65.2, 0.5e+3, 0.2, solid_fuel_req, 250e+3, 250e+3, 240, 225, [kerbalparts.solid_fuel(433)]))
kerbalparts.add(EnginePartType('rmbacc', 'Rockomax BACC Solid Fuel Booster', ['small', 'small'], 190, 1.5e+3, 0.2, solid_fuel_req, 315e+3, 315e+3, 250, 230, [kerbalparts.solid_fuel(850)]))
kerbalparts.add(EnginePartType('srb', 'S1 SRB-KD25K', ['small', 'small'], 300, 3e+3, 0.2, solid_fuel_req, 650e+3, 650e+3, 250, 230, [kerbalparts.solid_fuel(2500)]))

del liquid_fuel_req

del Catalogue, ResourceType, CrewPodPartType, ProbePodPartType, FuelTankPartType, EnginePartType
del FeedPartType

