from ._framednode import FramedNode
from .._frame import surface_frame, parse_surface_frame
from .._vector import statevector


class StaticSite(FramedNode):
    def __init__(self, parent_node, lat, lon, alt):
        FramedNode.__init__(self, parent_node, surface_frame(parent_node.surface, lat, lon, alt))
    
    def statevector(self, t):
        return self.frame.toinertial(statevector.zero(), t)

