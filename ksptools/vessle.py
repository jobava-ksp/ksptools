from . import part


class Stage(object):
    def __init__(self, parts, next=[]):
        self.next = next
        self.part = part.group(parts)
    
    def total_mass(self, state):
        return self.part.mass + sum(n.total_mass(state) for n in self.next)
    
    def total_dragcoef(self, state):
        return self.part.dragcoef + sum(n.total_mass(state) for n in self.next)
    
    def mass(self, state):
        raise NotImplementedError
    
    def dragcoef(self, state):
        raise NotImplementedError
