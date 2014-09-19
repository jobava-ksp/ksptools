from ._body import Body


class CelestialBody(Body):
    def __init__(self, parent_node):
        Body.__init__(parent_node)

