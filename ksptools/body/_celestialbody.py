from ._body import Body
from ._locallity import OrbitalFrame

class CelestialBody(Body):
    def __init__(self, parent_node, kepler):
        Body.__init__(parent_node, OrbitalFrame(kepler))

