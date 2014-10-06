import re
from numpy import array, cos, cosh, mat, pi, sin, sinh, sqrt, zeros
from numpy.linalg import norm


uniti = array([1,0,0])
unitj = array([0,1,0])
unitk = array([0,0,1])


def rotx(t):
    return mat([[1,      0,      0],
                [0, cos(t), -sin(t)],
                [0, sin(t),  cos(t)]])

def roty(t):
    return mat([[ cos(t), 0, sin(t)],
                [      0, 1,       0],
                [-sin(t), 0, cos(t)]])

def rotz(t):
    return mat([[ cos(t), -sin(t), 0],
                [ sin(t),  cos(t), 0],
                [      0,       0, 1]])
 
def rotaxis(axis,t):
    x,y,z = axis
    ct, st = cos(t), sin(t)
    return mat([[ct+x*x*(1-ct),   x*y*(1-ct)-z*st, x*z*(1-ct)+y*st],
                [y*x*(1-ct)+z*st, ct+y*y*(1-ct),   y*z*(1-ct)-x*st],
                [z*x*(1-ct)-y*st, z*y*(1-ct)+x*st, ct+z*z*(1-ct)  ]])

def rotzxz(a,b,c):
    return rotz(c)*rotx(b)*rotz(a)

def unit(v):
    if norm(v) == 0:
        return zeros(3)
    return v/norm(v)

def _limit(low, high):
    def decfunc(func):
        def limited_func(x):
            return func(min(high, max(low, x)))
        limited_func.__name__ = func.__name__
        limited_func.__doc__ = func.__doc__
        return limited_func
    return decfunc

@_limit(-5.0477e+5, 1.0e+200)
def C(z):
    if z > 0:
        return (1 - cos(z**0.5))/z
    elif z < 0:
        return (cosh(sqrt(-z)) - 1)/-z
    return 0.5

@_limit(-5.0477e+5, 1.0e+200)
def S(z):
    if z > 0:
        return (sqrt(z) - sin(z**0.5)) / sqrt(z)**3
    elif z < 0:
        return (sinh(sqrt(-z)) - sqrt(-z)) / sqrt(-z)**3
    return 1.0 / 6.0


class UnitMetaClass(type):
    def __getitem__(cls, key):
        return cls._units[key]

    def __setitem__(cls, key, value):
        cls._units[key] = value


class Unit(object):
    __metaclass__ = UnitMetaClass

    def __init__(self, name):
        self.name = name
        self._cfa = {name: 1}
        self._cfb = {name: 0}
    
    @classmethod
    def parse(cls, strvalue):
        unit_regex = '(?P<unit>({}))?'.format(')|('.join(v.name for v in cls._units.values()))
        value_regex = r'(?P<value>(\-|\+)?[0-9]+(\.[0-9]+)?(e(\+|\-)?[0-9]+)?)'
        m = re.match(value_regex + '\s*' + unit_regex, strvalue.strip())
        if not m:
            print(strvalue)
        if not m.group('unit'):
            return float(m.group('value')), None
        else:
            return float(m.group('value')), cls[m.group('unit')]
    
    @classmethod
    def parseto(cls, strvalue, unitname):
        value, srcunit = cls.parse(strvalue)
        if srcunit is None:
            return value
        else:
            destunit = cls[unitname]
            return cls.convert(value, srcunit, destunit)
    
    @classmethod
    def convert(cls, value, srcunit, destunit):
        return value * cls[srcunit.name]._cfa[destunit.name] + cls[srcunit.name]._cfb[destunit.name]
    
    @classmethod
    def addunit(cls, name, cfdict=None):
        if name not in cls._units:
            cls._units[name] = cls(name)
        if cfdict is not None:
            cls._addcfmap(name, cls[name], cfdict)
        return cls[name]
    
    @classmethod
    def _addcfmap(cls, name, unit, cfdict):
        for cfname, (cfa, cfb) in cfdict.items():
            cfunit = cls.addunit(cfname)
            cls.addcf(unit, cfunit, cfa, cfb)
        cfname_list = list(cfdict.keys())
        for bindex in range(len(cfname_list)-1):
            bname = cfname_list[bindex]
            bunit = cls[bname]
            for cname in cfname_list[bindex:]:
                cunit = cls[cname]
                Abc = unit._cfa[cname] * bunit._cfa[name]
                Bbc = unit._cfa[cname] * bunit._cfb[name] + unit._cfb[cname]
                cls.addcf(bunit, cunit, Abc, Bbc)
    
    @classmethod
    def addcf(cls, unit, cfunit, cfa, cfb):
        unit._cfa[cfunit.name], unit._cfb[cfunit.name] = cfa, cfb
        cfunit._cfa[unit.name], cfunit._cfb[unit.name] = 1/cfa, -cfb/cfa
    
    def __str__(self):
        return '{}'.format(self.name)
    
    def __repr__(self):
        return '<{}:{}>'.format(type(self).__name__, str(self))
    
    _units = dict()


Unit.addunit('m',
    {'pm': (1.0e+12, 0),
     'nm': (1.0e+9, 0),
     'um': (1.0e+6, 0),
     'mm': (1.0e+3, 0),
     'Km': (1.0e-3, 0),
     'Mm': (1.0e-6, 0),
     'Gm': (1.0e-9, 0),
     'Tm': (1.0e-12,0)})
Unit.addunit('deg',
    {'rad': (pi/180.0, 0)})
Unit.addunit('sec',
    {'psec': (1.0e+12, 0),
     'nsec': (1.0e+9, 0),
     'usec': (1.0e+6, 0),
     'msec': (1.0e+3, 0),
     'min': (1/60.0, 0),
     'hour': (1/3600.0, 0),
     'deg': (1/3600.0,0)})
Unit.addunit('K',
    {'C': (1.0, -273.15),
     'F': (1.8, -459.67),
     'R': (1.8, 0)})


def asunit(value, unit):
    return Unit.parseto(value, unit)


def asunits(values, units):
    return array([Unit.parseto(v,u) for v, u in zip(values, units)])


 
