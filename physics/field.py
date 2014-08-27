import scipy.constants as const

from numpy import dot
from numpy.linalg import norm

from .node import ModelNode

class Field(ModelNode):
    def __init__(self, modelset=['any']):
        ModelNode.__init__(self, modelset)
    
    def isbounded(self, body):
        raise NotImplementedError
    
    def force(self, body):
        raise NotImplementedError


class CenteredField(Field):
    def __init__(self, body, modelset=['any'], soi=None):
        Field.__init__(self, modelset)
        self.center = body
        self.soi = soi
    
    def isbounded(self, body):
        if self.soi is not None:
            return norm(body.x - center.x) < self.soi
        else:
            return True


class GravityField(CenteredField):
    def __init__(self, rigidbody, soi=None):
        CenteredField.__init__(rigidbody, ['newton', 'n-body'], soi)
        self.std_g_param = rigidbody.mass * const.G
    
    def force(self, body):
        r = body.x - self.center.x
        return (-(r/norm(r))*(self.std_g_param * body.mass)) / dot(r,r)


class AtmosphericField(CenteredField):
    def __init__(self, rigidbody, surface_p, surface_height, scale_height, soi):
        CenteredField.__init__(rigidbody, ['newton', 'n-body', 'kepler'], soi)
        self.surface_p = surface_p
        self.surface_height = surface_height
        self.scale_height = scale_height
    
    def force(self, body):
        r = body.x - self.center.x
        v = body.v - self.center.v
        alt = norm(r) - self.surface_height
        p = self.surface_p * const.e**(-alt/self.scale_height)
        d = v/norm(v)
        CdA = body.coefdrag(d) * body.area(d)
        return d * 0.5 * p * dot(v,v) * CdA


