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
        if not line.startswith("#"):
            x.append(float(line.split(',')[0]))
            y.append(float(line.split(',')[1]))
            u.append(float(line.split(',')[2]))
            v.append(float(line.split(',')[3]))
            P.append(float(line.split(',')[4]))
            omega.append(float(line.split(',')[5]))

    return np.array(x),np.array(y),np.array(u),np.array(v),np.array(P),np.array(omega)

filename = str(sys.argv[1])
filename_read = filename + ".dat"
velocity = str(sys.argv[2])
startX = int(sys.argv[3])
startY = int(sys.argv[4])
width = int(sys.argv[5])-2
height = int(sys.argv[6])-2
dx = float(sys.argv[7])
separation = float(sys.argv[8])+1

t1_minX = startX * dx
t1_maxX = dx * (startX + width) - t1_minX
t1_minY = startY * dx
t1_maxY = dx * (startY + height) - t1_minY

t2_minX = dx * (startX + width + separation)
t2_maxX = dx * (startX + 2 * width + separation) - t2_minX
t2_minY = t1_minY
t2_maxY = t1_maxY


x,y,u,v,P,omega = ParseFile(filename_read)

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

plotfile_vorticity_streamfunction = "vorticity_streamfunction_" + filename + ".pdf"
plotfile_velocityx = "x-velocity_" + filename + ".pdf"
plotfile_velocityy = "y-velocity_" + filename + ".pdf"
plotfile_pressure = "pressure_" + filename + ".pdf"


matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'

# Plots the vorticity and the streamfunction
title = "Vorticity and Streamfunction"
plt.clf()
plt.rc("text", usetex=True)
plt.rc("font", family="serif")
CS = plt.contourf(xi,yi,zi,200,cmap=plt.cm.jet)
CB = plt.colorbar(CS)
CB.set_label("vorticity $(\omega)$")
plt.xlabel("$x$")
plt.ylabel("$y$")
plt.title(title)
currentAxis = plt.gca()
currentAxis.add_patch(Rectangle((t1_minX,t1_minY),t1_maxX,t1_maxY,hatch='/',facecolor="white"))#,alpha=0.75))fill=False
#currentAxis.add_patch(Rectangle((t2_minX,t2_minY),t2_maxX,t2_maxY,hatch='/',facecolor="white"))#,alpha=0.75))fill=False
plt.streamplot(xm, ym, u2, v2, density=(1,1), color='k')#, linewidth=lw)
plt.ylim([0.3,1.5])
plt.xlim([0.3,1.5])
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
currentAxis = plt.gca()
currentAxis.add_patch(Rectangle((t1_minX,t1_minY),t1_maxX,t1_maxY,hatch='/',facecolor="white"))#,alpha=0.75))fill=False
#currentAxis.add_patch(Rectangle((t2_minX,t2_minY),t2_maxX,t2_maxY,hatch='/',facecolor="white"))#,alpha=0.75))fill=False
plt.xlabel("$x$")
plt.ylabel("$y$")
plt.title(title)
plt.ylim([0.3,1.5])
plt.xlim([0.3,1.5])
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
currentAxis = plt.gca()
currentAxis.add_patch(Rectangle((t1_minX,t1_minY),t1_maxX,t1_maxY,hatch='/',facecolor="white"))#,alpha=0.75))fill=False
#currentAxis.add_patch(Rectangle((t2_minX,t2_minY),t2_maxX,t2_maxY,hatch='/',facecolor="white"))#,alpha=0.75))fill=False
plt.xlabel("$x$")
plt.ylabel("$y$")
plt.title(title)
plt.ylim([0.3,1.5])
plt.xlim([0.3,1.5])
plt.savefig(plotfile_velocityy, format="pdf", bbox_inches="tight")

# Plots the pressure
title = "Pressure and Streamfunction"# $(R_{e} = " + Re + ")$"
zi = griddata(x, y, P, xi, yi, interp="linear")
plt.clf()
plt.rc("text", usetex=True)
plt.rc("font", family="serif")
CS = plt.contourf(xi,yi,zi,200,cmap=plt.cm.jet)
CB = plt.colorbar(CS)
CB.set_label("Pressure $(P)$")
plt.streamplot(xm, ym, u2, v2, density=(1,1), color='k')#, linewidth=lw)
currentAxis = plt.gca()
currentAxis.add_patch(Rectangle((t1_minX,t1_minY),t1_maxX,t1_maxY,hatch='/',facecolor="white"))#,alpha=0.75))fill=False
#currentAxis.add_patch(Rectangle((t2_minX,t2_minY),t2_maxX,t2_maxY,hatch='/',facecolor="white"))#,alpha=0.75))fill=False
plt.xlabel("$x$")
plt.ylabel("$y$")
plt.title(title)
plt.ylim([0.3,1.5])
plt.xlim([0.3,1.5])
plt.savefig(plotfile_pressure, format="pdf", bbox_inches="tight")
