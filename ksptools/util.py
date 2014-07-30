

def rotx(t):
    from numpy import mat, sin, cos
    return mat([[1,       0,      0],
                [0,  cos(t), sin(t)],
                [0, -sin(t), cos(t)]])

def roty(t):
    from numpy import mat, sin, cos
    return mat([[cos(t), 0, -sin(t)],
                [0,      1,       0],
                [sin(t), 0,  cos(t)]])

def rotz(t):
    from numpy import mat, sin, cos
    return mat([[cos(t), -sin(t), 0],
                [sin(t),  cos(t), 0],
                [     0,       0, 1]])

def veccos(a,b):
    from numpy import dot
    from numpy.linalg import norm
    if norm(a)*norm(b) == 0:
        return 0.0
    return dot(a,b)/(norm(a)*norm(b))

def vecsin(a,b,zisup=True):
    from numpy import array, cross
    from numpy.linalg import norm
    if len(a) == 2:
        return veccos(array([-a[1], a[0]]), b)
    elif len(a) == 3:
        v = cross(a,b)/(norm(a)*norm(b))
        if zisup and v[2] < 0:
            return -norm(v)
        else:
            return norm(v)

def arcvec(a,b):
    from numpy import arcsin, arccos, pi
    cost = veccos(a,b)
    sint = vecsin(a,b)
    if sint >= 0.0:
        return arccos(cost)
    else:
        return 2*pi - arccos(cost)

def Ax(A,x):
    from numpy import mat
    return (A*mat(x).T).A1

def project(a,b):
    from numpy import dot
    return b*(dot(a,b)/dot(b,b))

def reject(a,b):
    from numpy import dot
    return a - b*(dot(a,b)/dot(b,b))
    
def unit(a):
    from numpy.linalg import norm
    return a/norm(a)

def projmat(*basis):
    from numpy import mat
    return mat([unit(v) for v in basis])

### functions for orbital stuffs

def rofaet(a,e,t):
    from numpy import cos
    return a*(1-e**2)/(1+e*cos(t))

def aofret(r,e,t):
    from numpy import cos
    if e == 1.0:
        return r/2.
    return r*(1+e*cos(t))/(1-e**2)

def Eofet(e,t):
    from numpy import arccos, sin, cos, pi
    if e == 1.0:
        return t
    a = (e+cos(t))/(1+e*cos(t))
    a = min(1.0, max(-1.0, a))
    return arccos(a)

def Mofet(e,t):
    from numpy import sin, pi
    if t >= 2*pi:
        return Mofet(e,t-2*pi) + 2*pi
    E = Eofet(e,t)
    m = E - e*sin(E)
    if t > pi:
        return 2*pi-m
    else:
        return m

def nofretu(r,e,t,u):
    from numpy import sqrt
    return sqrt(u/aofret(r,e,t)**3)

def drofaet(a,e,t):
    from numpy import cos, sin
    da = (1-e**2)/(1+e*cos(t))
    de = a*((e**2+1)*cos(t) + 2*e)/((e*cos(t)+1)**2)
    dt = a*e*(e**2-1)*sin(t)/((e*cos(t)+1)**2)
    return da, de, dt

def daofret(r,e,t):
    from numpy import sin, cos
    if e == 1.:
        dr = float('-inf')
        de = float('-inf')
        dt = 0.0
    else:
        dr = (e+cos(t))/(1-e**2)
        de = r*((e**2+1)*cos(t)+2*e)/((e**2-1)**2)
        dt = e*r*sin(t)/(e**2-1)
    return dr, de, dt

def dEofet(e,t):
    from numpy import sqrt, sin, cos, pi
    
    if sin(t) == 0.:
        de = 0.0
        dt = 0.0
    elif e == 1.:
        de = 0.0
        dt = 0.0
    else:
        g = sqrt(-(e**2-1)*sin(t)**2/((e*cos(t)+1)**2))
        de = g/(e**2-1)
        dt = (1.0/sin(t))*g
    return de, dt

def dMofet(e,t):
    from numpy import cos, sin
    E = Eofet(e,t)
    dE_e, dE_t = dEofet(e,t)
    de = dE_e*(1-e*cos(E)) - sin(E)
    dt = dE_t*(1-e*cos(E))
    return de, dt

def dnofretu(r,e,t,u):
    from numpy import sqrt
    da_r, da_e, da_t = daofret(r,e,t)
    g = -(3./2.)*sqrt(u/aofret(r,e,t)**5)
    dr = g*da_r
    de = g*da_e
    dt = g*da_t
    du = 0.0
    return dr, de, dt, du


    
