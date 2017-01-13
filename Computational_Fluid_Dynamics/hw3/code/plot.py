#!/usr/bin/python

from __future__ import print_function, division
import math
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import sys
import matplotlib.cm as cm
from matplotlib.mlab import griddata
import scipy.interpolate

def ParseFile(filename):
    x = []
    y = []
    omega = []
    psi = []
    u = []
    v = []
    length = 3

    infile = open(filename, 'r')

    lines = [line.strip() for line in infile.readlines()]

    for line in lines:
        x.append(float(line.split(',')[0]))
        y.append(float(line.split(',')[1]))
        omega.append(float(line.split(',')[2]))
        psi.append(float(line.split(',')[3]))
        u.append(float(line.split(',')[4]))
        v.append(float(line.split(',')[5]))

    return np.array(x),np.array(y),np.array(omega),np.array(psi),np.array(u),np.array(v)

#title = str(sys.argv[1])
filename = "finished.dat"
x,y,omega,psi,u,v = ParseFile(filename)
xi = np.linspace(min(x), max(x))
yi = np.linspace(min(y), max(y))
zi = griddata(x, y, omega, xi, yi, interp="linear")
xsize = np.sqrt(len(x))
ysize = np.sqrt(len(y))
xm = np.linspace(min(x), max(x), xsize)
ym = np.linspace(min(y), max(y), ysize)
speed = []
lw = []

u2 = np.reshape(u,(xsize,ysize))
v2 = np.reshape(v,(xsize,ysize))

speed = np.sqrt(u2*u2 + v2*v2)
lw = 5*speed / speed.max()

np.array(speed)

plotfile = "plot_vorticity.pdf"
plotfile2 = "plot_vorticity_streamfunction.pdf"

matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'

plt.clf()
CS = plt.contourf(xi,yi,zi,200,cmap=plt.cm.jet)
plt.colorbar()
plt.xlabel("x")
plt.ylabel("y")
plt.title("Vorticity and Streamfunction")
plt.streamplot(xm, ym, u2, v2, density=(1,1), color='k')#, linewidth=lw)
plt.savefig(plotfile, format="pdf")
