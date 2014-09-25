from ._framednode import FramedNode


class Body(FramedNode):
    def __init__(self, parent_node, frame):
        FramedNode.__init__(self, parent_node, frame)

