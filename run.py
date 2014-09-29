from __future__ import division

import pprint

from kerbolsystem import kerbolsystem as ksys
from ksptools import *

import numpy as np
import numpy.linalg as la

import scipy.constants as spconst





class stagedcontroller(simcontroller):
    def __init__(self):
        pass
    
    def _staged_update(self, stage_name, depleted, max_twr, state):
        raise NotImplementedError
    
    def _initialize(self, payload, stages):
        self._stage = stage('payload', payload, 0, 0, 0, 0, None)
        self._state = None
        macc = paylaod
        for i in range(len(stages)):
            name, me, mf, Tmax, isp_0, isp_1 = stages[i]
            self._stage = (name, me + macc, mf, Tmax, isp_0, isp_1, self._stage)
            macc += me + mf
    
    def _update(self, m_rv, body_ijk_lla, pav, t):
        stage_name, me, _, Tmax, isp_0, isp_1, next_stage = self._stage
        if m <= me:
            _, coefd, Tdir, drop = self._staged_update(stage_name, True, Tmax / (m*spconst.g), (m_rv, body_ijk_lla, pav))
            f = 0
        else:
            f, coefd, Tdir, drop = self._staged_update(stage_name, False, Tmax / (m*spconst.g), (m_rv, body_ijk_lla, pav))
        if drop:
            self._stage = next_stage
        
        return (Tmax, isp_0, isp_1, f), coefd, Tdir
    
    @staticmethod
    def stage(name, me, mp, Tmax, isp_0, isp_1):
        return name, me, mf, Tmax, isp_0, isp_1


