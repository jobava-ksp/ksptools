from __future__ import division

import pprint

from kerbolsystem import kerbolsystem as ksys
from ksptools import *

import numpy as np
import numpy.linalg as la




from ksptools._frame import surface_frame, geodetic_frame

surface = surface_frame(geodetic_frame(3, 3, 0, 0, 0, np.pi/4.0), 0, 0, 1)

