import numpy as np
from mayavi.mlab import *
from colorsys import hsv_to_rgb
import sys

def color(theta, phi):
    while(phi < 0):
        phi += 2.0*np.pi
    h = phi / (2.0 * np.pi)
    s = 0.5
    v = 1.0 - 0.9999999*(theta / np.pi)
    print("h = {}, s = {}, v ={}".format(h,s,v))
    return hsv_to_rgb(h,s,v)


def plot_mag(file):
    data = np.load(file)

    x = data['m_x']
    y = data['m_y']
    z = data['m_z']

    phi   = data['m_phi']
    theta = data['m_theta']

    u = np.sin(theta) * np.cos(phi)
    v = np.sin(theta) * np.sin(phi)
    w = np.cos(theta)

    for i in range(x.shape[0]):
        r,g,b = color(theta[i], phi[i])
        print("R: {}, G: {}, B: {}".format(r,g,b))
        obj = quiver3d(x[i], y[i], z[i], u[i], v[i], w[i],
                line_width=3, color=(r,g,b), colormap='hsv', 
                scale_factor=0.3, mode='arrow',resolution=25)

plot_mag("/home/matthias/STB/output/dbg.npz")
