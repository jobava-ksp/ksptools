from ._body import Body


class RigidBody(Body):
    def __init__(self, parent_node, frame):
        Body.__init__(parent_node, frame)

