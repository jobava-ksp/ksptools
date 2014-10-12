from ksptools import *


import numpy as np
import numpy.linalg as la


stageA = partialstage('A', 0.55e+3, 50e+3, 390, 300)
stageB = partialstage('B', 1.55e+3, 200e+3, 370, 320)
sol = minimizefuel(1.0e+3, [stageA, stageB], [1.4, 1.7], 5000)
print(sol)


