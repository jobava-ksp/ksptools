from __future__ import division

import pprint

import ksptools as ksp
import ksptools.body as kspbody
import ksptools.orbit as kspo
import ksptools.unit as kspunit
import ksptools.util as kspu
import ksptools.kerbaltime as kspt
import ksptools.transfer as ksptransfer
import ksptools.sim as kspsim

import numpy as np
import numpy.linalg as npla

#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D

from kerbolsystem import kerbolsystem


class LaunchController(kspsim.Controller):
    def __init__(self, vessle, body, lon, lat, epoch):
        kspsim.Controller.__init__(self, body.geocentric(lon, lat, epoch)) 


class Vessle(kspbody.Body):
    pass


def launch(controller, rocket):
    kspsim.run(rocket, controller)

kerbin = kerbolsystem['kerbin']
ksc_lat, ksc_lon = kspunit.degloc_to_radloc(-np.array([0,6,9]),-np.array([74,34,31]))

ksc_loc = kerbin.geocentric(ksc_lon, ksc_lat, 0)
ref_loc = kerbin.geocentric(0,0,0)


