from ._node import Node


class Field(FramedNode):
    def __init__(self, parent_node):
        Node.__init__(parent_node, parent_node.surface_frame)


