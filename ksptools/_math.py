from numpy import array, cos, cosh, mat, sin, sinh, sqrt

uniti = array([1,0,0])
unitj = array([0,1,0])
unitk = array([0,0,1])


def rotx(t):
    return mat([[1,       0,      0],
                 [0,  cos(t), sin(t)],
                 [0, -sin(t), cos(t)]])

def roty(t):
    return mat([[cos(t), 0, -sin(t)],
                 [     0, 1,       0],
                 [sin(t), 0,  cos(t)]])

def rotz(t):
    return mat([[ cos(t), sin(t), 0],
                 [-sin(t), cos(t), 0],
                 [      0,      0, 1]])

def rotzxz(a,b,c):
    return rotz(c)*rotx(b)*rotz(a)
    

def C(z):
    if z > 0:
        return (1 - cos(z**0.5))/z
    elif z < 0:
        return (cosh(sqrt(-z)) - 1)/-z
    return 0.5


def S(z):
    if z > 0:
        return (sqrt(z) - sin(z**0.5)) / sqrt(z)**3
    elif z < 0:
        return (sinh(sqrt(-z)) - sqrt(-z)) / sqrt(-z)**3
    return 1.0 / 6.0


