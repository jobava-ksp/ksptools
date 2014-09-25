from __future__ import division

import pprint

from kerbolsystem import kerbolsystem as ksys
from ksptools import *

kerbin = ksys['kerbin']
orbit = kepler(0.001, 6.9e+5, 0, 0, 0, kerbin.GM, 0)
print(solve_hohmann_transfer(kerbin, orbit, 0, 7.2e+5).totaldv())


