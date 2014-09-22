from ksptools._math import *
from testutil import *

import numpy as np


class Test__Math(TestCase):
    def setUp(self):
        load_testdata(self, 'data._math_data')
    
    def _array_eq(self, x, y):
        self.assertTrue(all(x.A1==y.A1))

    def test_rotx(self):
        runtest(self, self.rotx_data, rotx, self._array_eq)
    
    def test_roty(self):
        runtest(self, self.roty_data, roty, self._array_eq)
    
    def test_rotz(self):
        runtest(self, self.rotz_data, rotz, self._array_eq)
    
    def test_rotzxz(self):
        runtest(self, self.rotzxz_data, rotzxz, self._array_eq)


