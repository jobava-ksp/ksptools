from ._framednode import FramedNode


class Field(FramedNode):
    def __init__(self, parent_node):
        FramedNode.__init__(self, parent_node, parent_node.surface_frame)


