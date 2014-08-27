from __future__ import division
from numpy import array, pi, matrix
import re


def degloc_to_radloc(lon, lat):
    DR = matrix([pi/180, pi/(180*60), pi/(180*60*60)]).T
    LL = matrix([lon,lat])
    return (LL*DR).T.A1


class UnitMetaClass(type):
    def __getitem__(cls, value):
        return cls._units[value]


class Unit(object):
    __metaclass__ = UnitMetaClass
    
    def __init__(self, name):
        self.name = name
        self._cfa = {name: 1}
        self._cfb = {name: 0}
    
    @classmethod
    def parse(cls, strvalue):
        unit_regex = '(?P<unit>({}))'.format(')|('.join(v.name for v in cls._units.values()))
        value_regex = r'(?P<value>(\-|\+)?[0-9]+(\.[0-9]+)?(e(\+|\-)?[0-9]+)?)'
        m = re.match(value_regex + '\s*' + unit_regex, strvalue)
        return float(m.group('value')), cls[m.group('unit')]
    
    @classmethod
    def parseto(cls, strvalue, unitname):
        value, srcunit = cls.parse(strvalue)
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

def tounit(value, unitname):
    return Unit.parseto(value, unitname)
    
