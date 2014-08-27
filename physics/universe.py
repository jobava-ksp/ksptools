from .node import TreeNode


class Universe(TreeNode):
    def __init__(self, universal_time=0):
        self.time = universal_time
        self.bodies = list()
        self.field = list()
    
    def addbody(self, body):
        self.bodies.append(body)
        self.addnode(body)
    
    def addfield(self, field):
        self.fields.append(field)
        self.addnode(field)
    
    def step(self, dt):
        for index in range(len(self.bodies)):
            body = self.bodies[index]
            ...
        
