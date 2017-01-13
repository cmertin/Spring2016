#!/usr/bin/python

from __future__ import print_function, division
import math
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import sys
import matplotlib.cm as cm
from matplotlib.mlab import griddata
from matplotlib.patches import Rectangle
import scipy.interpolate

def ParseFile(filename):
    x = []
    y = []
    u = []
    v = []
    P = []
    omega = []

    infile = open(filename, 'r')

    lines = [line.strip() for line in infile.readlines()]

    for line in lines:
        x.append(float(line.split(',')[0]))
        y.append(float(line.split(',')[1]))
        u.append(float(line.split(',')[2]))
        v.append(float(line.split(',')[3]))
        P.append(float(line.split(',')[4]))
        omega.append(float(line.split(',')[5]))

    return np.array(x),np.array(y),np.array(u),np.array(v),np.array(P),np.array(omega)

if(len(sys.argv) == 1):
    Re = str(100)
else:
    Re = str(str(sys.argv[1]))
#Re = str(sys.argv[1])
filename = "Re=" + Re + "_data.dat"
x,y,u,v,P,omega = ParseFile(filename)
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

plotfile_vorticity_streamfunction = "vorticity_streamfunction_" + Re + ".pdf"
plotfile_velocityx = "x-velocity_" + Re + ".pdf"
plotfile_velocityy = "y-velocity_" + Re + ".pdf"
plotfile_pressure = "pressure_" + Re + ".pdf"


matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'

# Plots the vorticity and the streamfunction
title = "Vorticity and Streamfunction"# $(R_{e} = " + Re + ")$"
plt.clf()
plt.rc("text", usetex=True)
plt.rc("font", family="serif")
CS = plt.contourf(xi,yi,zi,200,cmap=plt.cm.jet)
CB = plt.colorbar(CS)
CB.set_label("vorticity $(\omega)$")
plt.xlabel("$x$")
plt.ylabel("$y$")
plt.title(title)
#currentAxis = plt.gca()
#currentAxis.add_patch(Rectangle((0.2,0.2),.2,.2,hatch='/',facecolor="white",alpha=0.75))#fill=False
plt.streamplot(xm, ym, u2, v2, density=(1,1), color='k')#, linewidth=lw)
plt.savefig(plotfile_vorticity_streamfunction, format="pdf", bbox_inches="tight")

# Plots the x velocity with vector arrows
title = "Horizontal Velocity"# $(R_{e} = " + Re + ")$"
zi = griddata(x, y, u, xi, yi, interp="linear")
plt.clf()
plt.rc("text", usetex=True)
plt.rc("font", family="serif")
CS = plt.contourf(xi,yi,zi,200,cmap=plt.cm.jet)
CB = plt.colorbar(CS)
CB.set_label("velocity $(u)$")
plt.quiver(x,y,u2, v2, color="white")
plt.xlabel("$x$")
plt.ylabel("$y$")
plt.title(title)
plt.savefig(plotfile_velocityx, format="pdf", bbox_inches="tight")

# Plots the y velocity with vector arrows
title = "Vertical Velocity"# $(R_{e} = " + Re + ")$"
zi = griddata(x, y, v, xi, yi, interp="linear")
plt.clf()
plt.rc("text", usetex=True)
plt.rc("font", family="serif")
CS = plt.contourf(xi,yi,zi,200,cmap=plt.cm.jet)
CB = plt.colorbar(CS)
CB.set_label("velocity $(v)$")
plt.quiver(x,y,u2, v2, color="white")
plt.xlabel("$x$")
plt.ylabel("$y$")
plt.title(title)
plt.savefig(plotfile_velocityy, format="pdf", bbox_inches="tight")

# Plots the pressure
title = "Pressure"# $(R_{e} = " + Re + ")$"
zi = griddata(x, y, P, xi, yi, interp="linear")
plt.clf()
plt.rc("text", usetex=True)
plt.rc("font", family="serif")
CS = plt.contourf(xi,yi,zi,200,cmap=plt.cm.jet)
CB = plt.colorbar(CS)
CB.set_label("Pressure $(P)$")
plt.xlabel("$x$")
plt.ylabel("$y$")
plt.title(title)
plt.savefig(plotfile_pressure, format="pdf", bbox_inches="tight")
