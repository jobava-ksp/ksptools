

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

def vecsin(a,b):
    from numpy import array, cross
    from numpy.linalg import norm
    if len(a) == 2:
        return veccos(array([-a[1], a[0]]), b)
    elif len(a) == 3:
        return norm(corss(a,b))/(norm(a)*norm(b))

def arcvec(a,b):
    from numpy import arcsin, arccos, pi
    cost = veccos(a,b)
    sint = vecsin(a,b)
    if sint >= 0:
        return arccos(cost)
    else:
        return 2*pi - arccos(cost)

def Ax(A,x):
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
