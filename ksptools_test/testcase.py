import unittest
from numpy.testing import *


class DataTestCase(unittest.TestCase):
    def setUp(self):
        from os.path import join as joinpath
        gdict = globals()
        execfile(joinpath('./ksptools_test/data/', type(self).dataname()), gdict)
        self._testsets = dict()
        self._testobj = dict()
        
        for k, v in [(l,w) for l,w in gdict.items() if l.startswith('testset_')]:
            self._testsets[k[8:]] = v
        for k, v in [(l,w) for l,w in gdict.items() if l.startswith('obj_')]:
            self._testobj[k[4:]] = v
    
    def run_testset(self, testset_name, testfunc, predicate):
        testset = self._testsets[testset_name]
        for exp_val, res_val in iter((y, testfunc(*x)) for x, y in testset):
            predicate(exp_val, res_val)
    
    def __getitem__(self, keyname):
        return self._testobj[keyname]
    
    @classmethod
    def dataname(cls):
        from os.path import basename
        return basename(__file__).split('.')[0][4:] + '_data'


