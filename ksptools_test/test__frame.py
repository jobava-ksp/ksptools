from ksptools_test.testcase import *


class TestFrame(DataTestCase):
    def test_surface_ijk(self):
        self.run_testset('surface_ijk', self['surface_ijk_func'], assert_almost_equal)
    
    def test_surface_inertial_statevector(self):
        self.run_testset('surface_inertial_statevector', self['surface_inertial_statevector_func'], assert_almost_equal)
    
    def test_geodetic_llav(self):
        self.run_testset('geodetic_llav', self['geodetic_llav_func'], assert_almost_equal)
    
    @classmethod
    def dataname(cls):
        return '_frame_data.py'


