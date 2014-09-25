from numpy.linalg import norm


class State(object):
    def __init__(self, ref, stv, t):
        self.vector = stv
        self.reference = ref
        self.time = t


class PathNode(object):
    def __init__(self, refi, stvi, ti, reff, stvf, tf):
        self.state_i = State(refi, stvi, ti)
        self.state_f = State(reff, stvf, tf)
    
    def totaldv(self):
        raise NotImplementedError


class Coast(PathNode):
    def __init__(self, ref, kepler, ti, tf=None):
        if tf is not None:
            PathNode.__init__(self,
                    ref, kepler.statevector_by_time(ti), ti,
                    ref, kepler.statevector_by_time(tf), tf)
        else:
            PathNode.__init__(self,
                    ref, kepler.statevector_by_time(ti), ti,
                    ref, None, None)
        
        self.orbit = kepler
    
    def totaldv(self):
        return 0


class Impulse(PathNode):
    def __init__(self, ref, stvi, stvf, ti):
        PathNode.__init__(self,
            ref, stvi, ti,
            ref, stvf, ti)
    
    def _get_dv(self):
        return self.state_f.vector.v - self.state_i.vector.v
    
    def totaldv(self):
        return norm(self._get_dv())
    
    dv = property(_get_dv)


class Path(PathNode):
    def __init__(self, path_list):
        si = path_list[0].state_i
        sf = path_list[-1].state_f
        PathNode.__init__(self,
            si.reference, si.vector, si.time,
            sf.reference, sf.vector, sf.time)
        self.path_list = path_list

    @staticmethod
    def coast(ref, kepler, ti, tf=None):
        return Coast(ref, kepler, ti, tf)
    
    @staticmethod
    def impulse(ref, stvi, stvf, t):
        return Impulse(ref, stvi, stvf, t)
    
    def totaldv(self):
        return sum([p.totaldv() for p in self.path_list])
    
