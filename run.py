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
stages = [stage('A', 0.675e+3, 1.0e+3, 50e+5, 390, 300)]

controller = LaunchController()
print(controller.sim(0.94, stages, ksys['ksc'].statevector(0), ksys['kerbin'], 0, 3*60))

