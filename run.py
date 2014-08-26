from __future__ import division

import pprint

import ksptools as ksp
import ksptools.body as kspbody
import ksptools.orbit as kspo
import ksptools.util as kspu
import ksptools.kerbaltime as kspt
import ksptools.transfer as ksptransfer
import ksptools.sim as kspsim

import numpy as np
import numpy.linalg as npla

#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D

from kerbolsystem import *


class LaunchController(kspsim.Controller):
    pass


class Vessle(kspbody.Body):
    pass


def launch(controller, rocket):
    kspsim.run(rocket, controller)