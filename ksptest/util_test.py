import unittest
import itertools
import numpy as np
import numpy.random as nprandom
import numpy.testing as nptest

from ksptools.util import *
from numpy import linspace, array, sin, cos, pi
from numpy.linalg import norm


class UtilTestCase(unittest.TestCase):
    def setUp(self):
        self.thetarange = linspace(0, 8*pi, 64, False)
        self.sinrange = sin(self.thetarange)
        self.cosrange = cos(self.thetarange)
        self.thetacossin = zip(self.thetarange, self.cosrange, self.sinrange)
        self.vectors = array([nprandom.rand(3) for i in range(64)])
    
    @unittest.skip("Not tested")
    def test_Ax(self):
        pass
    
    def iter_vector_theta(self, func):
        for t, ct, st in self.thetacossin:
            for a in self.vectors:
                yield a, func(a, t), t, ct, st
    
    def test_rotx(self):
        for a, b, t, ct, st in self.iter_vector_theta(lambda v, t: Ax(rotx(t),v)):
            ax, ay, az = a
            bx, by, bz = b
            nptest.assert_allclose(array([bx, by, bz]), array([ax, ay*ct - az*st, ay*st + az*ct]))
    
    def test_roty(self):
        for a, b, t, ct, st in self.iter_vector_theta(lambda v, t: Ax(roty(t),v)):
            ax, ay, az = a
            bx, by, bz = b
            nptest.assert_allclose(array([bx, by, bz]), array([ct*ax + st*az, ay, -st*ax + ct*az]))
    
    def test_roty(self):
        for a, b, t, ct, st in self.iter_vector_theta(lambda v, t: Ax(rotz(t),v)):
            ax, ay, az = a
            bx, by, bz = b
            nptest.assert_allclose(array([bx, by, bz]), array([ct*ax - st*ay, st*ax + ct*ay, az]))
    
    @unittest.skip("Not tested")
    def test_cosssin(self):
        pass
    
    def test_veccos(self):
        for t, ct, st in self.thetacossin:
            v = Ax(rotz(t), uniti)
            nct = veccos(uniti, v)
            nptest.assert_allclose(array([ct]), array([nct]))
    
    def test_vecsin(self):
        for t, ct, st in self.thetacossin:
            v = Ax(rotz(t), uniti)
            nst = vecsin(uniti, v)
            nptest.assert_allclose(array([st]), array([nst]))
    
    def test_arcvec(self):
        for t, ct, st in self.thetacossin:
            t = t % (2*pi)
            v = Ax(rotz(t), uniti)
            nt = arcvec(uniti, v)
            nptest.assert_allclose(array([nt]), array([t]))
    
    @unittest.skip("Not tested")
    def test_project(self):
        for v in self.vectors:
            vi = project(v, uniti)
            vj = project(v, unitj)
            vk = project(v, unitk)
            nptest.assert_allclose(array([norm(vi), abs(veccos(v,vi))]), array([abs(v[0]), 1.]))
            nptest.assert_allclose(array([norm(vj), abs(veccos(v,vj))]), array([abs(v[1]), 1.]))
            nptest.assert_allclose(array([norm(vk), abs(veccos(v,vk))]), array([abs(v[2]), 1.]))
    
    @unittest.skip("Not tested")
    def test_reject(self):
        for v in self.vectors:
            vi = reject(v, uniti)
            vj = reject(v, unitj)
            vk = reject(v, unitk)
            nptest.assert_allclose(array([norm(vi), abs(veccos(v,vi))]), array([nrom(v[1:]), 0.]))
            nptest.assert_allclose(array([norm(vj), abs(veccos(v,vj))]), array([norm(array([v[0], v[2]])), 0.]))
            nptest.assert_allclose(array([norm(vk), abs(veccos(v,vk))]), array([norm(v[:2]), 0.]))
    
    def test_unit(self):
        for v in self.vectors:
            u = unit(v)
            nptest.assert_allclose(array([norm(u), veccos(u,v)]), array([1.,1.]))
    
