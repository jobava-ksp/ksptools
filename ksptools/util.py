from __future__ import division
from numpy import array, sin, cos, mat, dot, cross, arccos, arcsin, arctan2, pi
from numpy.linalg import norm


uniti = array([1.,0.,0.])
unitj = array([0.,1.,0.])
unitk = array([0.,0.,1.])

def rotx(t):
    return mat([[1,       0,      0],
                 [0,  cos(t), sin(t)],
                 [0, -sin(t), cos(t)]])

def roty(t):
    return mat([[cos(t), 0, -sin(t)],
                 [0,      1,       0],
                 [sin(t), 0,  cos(t)]])

def rotz(t):
    return mat([[cos(t), -sin(t), 0],
                [sin(t),  cos(t), 0],
                [     0,       0, 1]])

def rotvec(u,t):
    x,y,z = u
    ct, st, _ = cossin(t)
    return mat([[ct+x*x*(1-ct),   x*y*(1-ct)-z*st, x*z*(1-ct)+y*st],
                [y*x*(1-ct)+z*st, ct+y*y*(1-ct),   y*z*(1-ct)-x*st],
                [z*x*(1-ct)-y*st, z*y*(1-ct)+x*st, ct+z*z*(1-ct)  ]])

def cossin(t, dim=3):
    return array([cos(t), sin(t)] + [0.]*(dim-2))

def veccos(a,b):
    _a = dot(a,b)/(norm(a)*norm(b))
    #if norm(a)*norm(b) == 0:
    #    return 0.0
    return min(1.0, max(-1.0, _a))

def vecsin(a,b,up=unitk):
    if len(a) == 2:
        return veccos(array([-a[1], a[0]]), b)
    elif len(a) == 3:
        v = cross(a,b)/(norm(a)*norm(b))
        if dot(v,up) < 0:
            return -norm(v)
        else:
            return norm(v)

def arcvec(a,b):
    cost = veccos(a,b)
    sint = vecsin(a,b)
    t = arctan2(sint, cost)
    if t < 0:
        t = 2*pi + t
    return t

def Ax(A,x):
    return (A*mat(x).T).A1

def project(a,b):
    return b*(dot(a,b)/dot(b,b))

def reject(a,b):
    return a - b*(dot(a,b)/dot(b,b))
    
def unit(a):
    return a/norm(a)

def projmat(*basis):
    return mat([unit(v) for v in basis])


class EulerAngle(object):    
    def __init__(self, p, t, s):
        Ap = rotz(p)  
        At = rotx(t)
        As = rotz(s)
        
        self.phi, self.theta, self.psi = p, t, s
        self.Ap, self.At, self.As = Ap, At, As
        self.A = Ap*At*As
    
    @classmethod
    def from_pts(cls,p,t,s):
        return cls(p,t,s)
    
    def rotate(self, r):
        return Ax(self.A,r)
    
    def derotate(self, r):
        return Ax(self.A.T,r)
    
    def __mul__(self, r):
        return Ax(self.A,r)
    
    def __div__(self, r):
        return Ax(self.A.T,r)


### -- matplotlib -- ###
def plotfunc(f, x0, x1, ax=None):
    from numpy import linspace
    import matplotlib.pyplot as plt
    x = linspace(x0, x1, 600)
    y = [f(t) for t in x]
    
    if ax is None:
        plt.plot(x,y)
        plt.show()
    else:
        ax.plot(x,y)


def plot_semi_orbit(kep, t0, t1, ax=None):
    from numpy import linspace
    import matplotlib.pyplot as plt
    
    t = linspace(t0,t1,600)
    r = [kep.r(i) for i in t]
    x, y, z = zip(*r)
    if ax is None:
        plotter = plt
    else:
        plotter = ax
    plotter.plot(x,y)

def plot_semi_orbit3d(kep, t0, t1, ax):
    from numpy import linspace
    import matplotlib.pyplot as plt
    
    t = linspace(t0,t1,600)
    r = [kep.r(i) for i in t]
    x, y, z = zip(*r)
    ax.plot(x,y,z)

def plot_rv(r, v, ax=None, scale=1.):
    from numpy import linspace
    import matplotlib.pyplot as plt
    
    if ax is None:
        plotter = plt
    else:
        plotter = ax
    rv = (r, r+scale*v)
    x, y, z = zip(*rv)
    plotter.plot(x,y)

def plot_rv3d(r, v, ax, scale=1.):
    from numpy import linspace
    import matplotlib.pyplot as plt

    rv = (r, r+v)
    x, y, z = zip(*rv)
    ax.plot(x,y,z)
    
