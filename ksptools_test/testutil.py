from unittest import TestCase

class DataTransformTest(object):
    def __init__(self, data_in, data_out):
        self.data_in = data_in
        self.data_out = data_out
    
    def run(self, func, predicate):
        predicate(self.data_out, func(*(self.data_in)))


def tcv(in_value, out_value):
    return DataTransformTest(in_value, out_value)


def runtest(tc, tcvlist, func, predicate=None):
    if predicate is None:
        predicate = tc.assertEqual
    for tcv in tcvlist:
        tcv.run(func, predicate)


def load_testdata(obj, datamodule):
    test_globals = dict()
    exec('from {} import *'.format(datamodule), test_globals)
    for k, v in test_globals.items():
        if not k.startswith('_'):
            setattr(obj, k, v)

