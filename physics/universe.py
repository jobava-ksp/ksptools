import itertools

from numpy import array
from .node import TreeNode
from . import force as forces

class Universe(TreeNode):
    def __init__(self, universal_time=0):
        self.time = universal_time
        self.bodies = list()
        self.fields = list()
        self.models = dict()
    
    def addbody(self, body, model):
        self.bodies.append(body)
        self.addnode(body)
        self.models[body] = model
    
    def addfield(self, field):
        self.fields.append(field)
        self.addnode(field)
    
    def setmodel(self, body, model):
        self.models[body] = model
    
    def prestep(self):
        self._forces = itertools.defaultdict(list)
    
    def applyforce_local(self, body, force, local_position=array([0,0,0])):
        self._forces[body].append(forces.force_local(body, force, local_position))
    
    def applyforce_global(self, body, force, global_position=array([0,0,0])):
        self._forces[body].append(forces.force_global(body, force, global_position))
    
    def applytorque_local(self, body, torque):
        self._forces[body].append(forces.torque_local(body, torque))
    
    def applytorque_global(self, body, torque):
        self._forces[body].append(forces.torque_global(body, torque))
    
    def _step_fields(self, dt):
        for index in range(len(self.bodies)):
            body = self.bodies[index]
            model = self.models[body]
            force = array([0,0,0])
            torque = array([0,0,0])
            for field in self.fields:
                if field.inmodel(model) and field.inbounds(body):
                    force += field.force(body)
            self.applyforce_local(body, force)
    
    def step(self, dt):
        self._step_fields(dt)
        for index in range(len(self.bodies)):
            body = self.bodies[index]
            force = sum(f for f in self._forces[body].force)
            torque = sum(t for t in self._forces[body].torque)
            model = self.models[body]
            if model in ['newton', 'n-body']:
                accel = force / body.mass
                accel_ang = torque / body.inertia(torque/norm(torque))
                body.transform_newton(dt, accel, accel_ang)
        self.time += dt
    
    def poststep(self):
        self._forces = None


