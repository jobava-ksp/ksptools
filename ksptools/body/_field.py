from ._framednode import FramedNode


class Field(FramedNode):
    def __init__(self, parent_node, frame):
        FramedNode.__init__(self, parent_node, frame)


