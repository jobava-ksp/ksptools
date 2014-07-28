import ksptools
import numpy as np
import numpy.linalg as la
import math
import matplotlib.pyplot as plt

from pprint import pprint


## name, radius, u, day, soi, atm(...), orbit(a,e,i,l,p,M,epoch)

qpi = 0.25*np.pi

sun = ksptools.body.CelestialBody('Sun', 1.0, 5.0, 0.0, 0.0, None, None)
a = ksptools.body.CelestialBody('a', 0.2, 0.04, 0.0, 0.0, None, ksptools.orbit.KeplerOrbit(sun, 6.0, 0.6, 0.0, 0.0, 0.0, 0.0, 0.0))
b = ksptools.body.CelestialBody('b', 0.2, 0.04, 0.0, 0.0, None, ksptools.orbit.KeplerOrbit(sun, 6.0, 0.6, qpi, 0.0, 0.0, 0.0, 0.0))
c = ksptools.body.CelestialBody('c', 0.2, 0.04, 0.0, 0.0, None, ksptools.orbit.KeplerOrbit(sun, 6.0, 0.6, 0.0, qpi, 0.0, 0.0, 0.0))
d = ksptools.body.CelestialBody('d', 0.2, 0.04, 0.0, 0.0, None, ksptools.orbit.KeplerOrbit(sun, 6.0, 0.6, qpi, 0.0, qpi, 0.0, 0.0))

def plot_orbit_xy(body):
    t, rx, ry, _, _, _, _ = get_orbit(body)
    plt.plot(rx,ry)

def plot_orbit_xz(body):
    t, rx, _, rz, _, _, _ = get_orbit(body)
    plt.plot(rx,rz)

def plot_orbit_xyz(ax, body):
    t, rx, ry, rz, _, _, _ = get_orbit(body)
    ax.plot(rx,ry,rz)

def get_orbit(body):
    T = 2*np.pi*np.sqrt(body.orbit.a**3/body.orbit.body.std_g_param)
    t = np.linspace(0,T,600,False)
    r, v = zip(*[body.orbit.rv(i) for i in t])
    rx, ry, rz = zip(*r)
    vx, vy, vz = zip(*v)
    return t, rx, ry, rz, vx, vy, vz

if __name__ == '__main__':
    from mpl_toolkits.mplot3d import Axes3D
    fig = plt.figure()
    fig.add_subplot(2, 2, 1)
    plot_orbit_xy(a)
    plot_orbit_xy(b)
    plot_orbit_xy(c)
    plot_orbit_xy(d)
    plt.axis('equal')
    fig.add_subplot(2, 2, 3)
    plot_orbit_xz(a)
    plot_orbit_xz(b)
    plot_orbit_xz(c)
    plot_orbit_xz(d)
    plt.axis('equal')
    ax = fig.add_subplot(2, 2, 2, projection='3d')
    plot_orbit_xyz(ax, a)
    plot_orbit_xyz(ax, b)
    plot_orbit_xyz(ax, c)
    plot_orbit_xyz(ax, d)
    plt.show()



