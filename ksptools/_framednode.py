from ._node import Node


class FramedNode(Node):
    def __init__(self, parent_node):
        Node.__init__(self, parent_node)

