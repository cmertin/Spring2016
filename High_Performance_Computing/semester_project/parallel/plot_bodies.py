#!/usr/bin/python

from __future__ import division, print_function
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import sys

def ReadData(filename):
    data = []
    mass = []
    x = []
    y = []
    z = []

    infile = open(filename, 'r')

    lines = [line.strip() for line in infile.readlines()]

    lines = lines[2:]

    for line in lines:
        mass.append(float(line.split(',')[0]))
        x.append(float(line.split(',')[1]))
        y.append(float(line.split(',')[2]))
        z.append(float(line.split(',')[3]))

    return mass,x,y,z

if len(sys.argv) == 1:
    size = str(1000)
else:
    size = str(sys.argv[1])

filename = size + "_bodies.dat"

mass,x,y,z = ReadData(filename)

cmhot = plt.cm.get_cmap("hot")
fig = plt.figure()
ax = fig.add_subplot(111,projection="3d")
l = ax.scatter(x,y,z,c=mass,cmap=cmhot)#,cmap=cmap)#,norm=norm)
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")
ax.set_xlim([-10,10])
ax.set_ylim([-10,10])
ax.set_zlim([-10,10])
ax.set_title("Celestial Bodies")
CB = fig.colorbar(l)
CB.set_label("Mass (Solar Masses)")
plt.show()
