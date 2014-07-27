import math
import re

_cf = dict({
        'rad': dict({
            'rad': 1,
            'deg': 180.0/math.pi
        }),
        'deg': dict({
            'deg': 1,
            'rad': math.pi/180.0
        }),
        'C': dict({
            'C': 1
        })
    })

def fromtounit(f, srcunit, destunit):
    return _cf[srcunit][destunit]*f

def tounit(expr, destunit):
    srcunit = re.search(r'(?P<unit>[^0-9]+)$',expr).group('unit')
    return fromtounit(float(expr[0:-len(srcunit)]), srcunit, destunit)

