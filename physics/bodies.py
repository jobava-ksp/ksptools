from __future__ import division
from . import body


class RigidKSPAprox(body.RigidBody):
    def __init__(self, position, velocity, mass, dragcoef):
        body.RigidBody.__init__(self, position, velocity, mass)
        self._cf = dragcoef
    
    def _coefdrag_local(self, vdir):
        return self._cf
    
    def _area_local(self, vdir):
        return 8e-3*self.mass
    
    def _inertia_local(self, axis):
        return (2/5)*(self.mass**2)*8e-3

