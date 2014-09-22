from __future__ import division
from testutil import *
import itertools
import numpy as np
import numpy.random as nprand

_range = 64
_theta = nprand.choice([(i/_range)*np.pi for i in range(_range)], 5)
_three_theta = list(itertools.product(_theta, _theta, _theta))

_cos_theta = [np.cos(t) for t in _theta]
_sin_theta = [np.sin(t) for t in _theta]

def _rotx(ct, st):
    return np.mat([[1,0,0],[0,ct,st],[0,-st,ct]])

def _roty(ct, st):
    return np.mat([[ct,0,-st],[0,1,0],[st,0,ct]])

def _rotz(ct, st):
    return np.mat([[ct,st,0],[-st,ct,0],[0,0,1]])

rotx_data = [tcv((t,), _rotx(ct,st)) for t, ct, st in zip(_theta, _cos_theta, _sin_theta)]
roty_data = [tcv((t,), _roty(ct,st)) for t, ct, st in zip(_theta, _cos_theta, _sin_theta)]
rotz_data = [tcv((t,), _rotz(ct,st)) for t, ct, st in zip(_theta, _cos_theta, _sin_theta)]
rotzxz_data = [tcv((t), _rotz(np.cos(t[2]), np.sin(t[2])) * _rotx(np.cos(t[1]), np.sin(t[1])) * _rotz(np.cos(t[0]), np.sin(t[0]))) for t in _three_theta]
 
