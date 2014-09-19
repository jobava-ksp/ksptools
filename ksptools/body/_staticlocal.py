from ._framednode import FramedNode


class StaticSite(FramedNode):
    def __init__(self, parent_node):
        FramedNode.__init__(parent_node, parent_node.surface_frame)
        
