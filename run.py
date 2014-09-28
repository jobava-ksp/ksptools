from __future__ import division

import pprint

from kerbolsystem import kerbolsystem as ksys
from ksptools import *

import numpy as np
import numpy.linalg as la

kerbol = ksys['kerbol']
kerbin = ksys['kerbin']
mun = ksys['mun']
duna = ksys['duna']

stv = mun.tolocal(mun.statevector(0), 0)

def print_stv(stv):
    print("<{}*{},{}*{}>".format(la.norm(stv.r), stv.r/la.norm(stv.r), la.norm(stv.v), stv.v/la.norm(stv.v)))

print_stv(mun.relative_statevector(duna, stv, 0))

