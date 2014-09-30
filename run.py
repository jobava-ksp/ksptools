from __future__ import division

import pprint

from kerbolsystem import kerbolsystem as ksys
from ksptools import *

import numpy as np
import numpy.linalg as la


class LaunchController(StagedController):
    def __init__(self):
        StagedController.__init__(self)
    
    def _staged_update(self, stage_name, depleted, max_twr, state):
        m_rv, body_ijk_llav, pav, t = state
        k = body_ijk_llav[3]
        return 1, 0.2, k, depleted


stage = StagedController.stage
stages = [stage('A', 0.675e+3, 1.0e+3, 50e+3, 390, 300)]

controller = LaunchController()

stv0 = ksys['ksc'].statevector(0)
print(vars(ksys['ksc'].frame))
print(ksys['kerbin'].llav(stv0, 0))
print(str(stv0))
p, a, v = ksys['kerbin'].atmstate_by_statevector(stv0, 0)
print(str(statevector(stv0.r, v)))

#print(vars(ksys['ksc'].frame))
#print(ksys['kerbin'].llav(stv0, 0))
print(map(str, controller.sim(0.94e+3, stages, stv0, ksys['kerbin'], 0, 60)))

