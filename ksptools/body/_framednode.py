from ._node import Node


class FramedNode(Node):
    def __init__(self, parent_node, frame):
        Node.__init__(self, parent_node)
        self.frame = frame
    
    def tolocal(self, rv, t):
        return self.frame.tolocal(rv, t)
    
    def toinertial(self, rv, t):
        return self.frame.toinertial(rv, t)
    
    def relative_statevector(self, node, stv, t):
        return self.walk(
            node,
            stv,
            lambda p, next, x: p.toinertial(x, t),
            lambda p, next, x: p.tolocal(x, t))

