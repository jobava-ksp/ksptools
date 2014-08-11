

def runsimulation(controller, state, start_time, dt=1.0e-6, end_time=None):
    world_time = start_time
    while not controller.terminal(state, world_time, dt):
        controller.step(state, world_time, dt)
        world_time += dt
    return state

def Controller(object):
    def __init__(self):
        pass
    
    def terminal(self, state, world_time, dt):
        raise NotImplementedError
    
    def step(self, state, world_time, dt):
        raise NotImplementedError


