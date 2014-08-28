import itertools
import collections

from numpy import array
from numpy.linalg import norm
from .node import TreeNode
from . import force as forces


class Universe(TreeNode):
    def __init__(self, universal_time=0):
        """
        :type universal_time: float
        """
        TreeNode.__init__(self)
        self.time = universal_time
        self.bodies = list()
        self.fields = list()
        self.models = dict()
        self._forces = None
    
    def addbody(self, body, model):
        """
        :type body: physics.body.RigidBody
        :type model: str
        """
        self.bodies.append(body)
        self.addnode(body)
        self.models[body] = model
    
    def addfield(self, field):
        """
        :type field: physics.field.Field
        """
        self.fields.append(field)
        self.addnode(field)
    
    def setmodel(self, body, model):
        """
        :type body: physics.body.RigidBody
        :type model: str
        """
        self.models[body] = model
    
    def prestep(self):
        self._forces = collections.defaultdict(list)
    
    def applyforce_local(self, body, force, local_position=array([0, 0, 0])):
        """
        :type body: physics.body.RigidBody
        :type force: numpy.ndarray
        :type local_position: numpy.ndarray
        """
        self._forces[body].append(forces.force_local(body, force, local_position))
    
    def applyforce_global(self, body, force, global_position=array([0, 0, 0])):
        """
        :type body: physics.body.RigidBody
        :type force: numpy.ndarray
        :type global_position: numpy.ndarray
        """
        self._forces[body].append(forces.force_global(force, global_position))
    
    def applytorque_local(self, body, torque):
        """
        :type body: physics.body.RigidBody
        :type torque: numpy.ndarray
        """
        self._forces[body].append(forces.torque_local(torque))
    
    def applytorque_global(self, body, torque):
        """
        :type body: physics.body.RigidBody
        :type torque: numpy.ndarray
        """
        self._forces[body].append(forces.torque_global(body, torque))
    
    def _step_fields(self):
        for index in range(len(self.bodies)):
            body = self.bodies[index]
            model = self.models[body]
            force = array([0, 0, 0])
            torque = array([0, 0, 0])
            for field in self.fields:
                if field.inmodel(model) and field.inbounds(body):
                    force += field.force(body)
            self.applyforce_local(body, force)
    
    def step(self, dt):
        """
        :type dt: float
        """
        self._step_fields()
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


