#!/usr/bin/python

from __future__ import division, print_function
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

def ReadData(filename):
    data = []
    x = []
    y = []
    z = []

    infile = open(filename, 'r')

    lines = [line.strip() for line in infile.readlines()]

    for line in lines:
        x.append(float(line.split(',')[1]))
        y.append(float(line.split(',')[2]))
        z.append(float(line.split(',')[3]))

    return x,y,z

sunfile = "sun.dat"
earthfile = "earth.dat"
marsfile = "mars.dat"

sunx,suny,sunz = ReadData(sunfile)
earthx,earthy,earthz = ReadData(earthfile)
marsx,marsy,marsz = ReadData(marsfile)


fig = plt.figure()
ax = fig.gca(projection="3d")
ax.plot(sunx,suny,sunz,label="Sun")
ax.plot(earthx,earthy,earthz,label="Earth")
ax.plot(marsx,marsy,marsz,label="Earth")
ax.legend()
#axplt.rc("text", usetex=True)
#plt.rc("font", family="serif")
ax.set_xlabel("$x$")
ax.set_ylabel("$y$")
ax.set_zlabel("$z$")
ax.set_title("Celestial Bodies")

'''
plt.clf()
plt.rc("text", usetex=True)
plt.rc("font", family="serif")
plt.xlabel("$x$")
plt.ylabel("$y$")
plt.zlabel("$z$")
plt.title("Celestial Bodies")
plt.plot(sunx,suny,sunz,label="Sun")
plt.plot(earthx,earthy,sunz,label="Earth")
plt.plot(marsx,marsy,sunz,label="Mars")
'''
plt.show()
    
