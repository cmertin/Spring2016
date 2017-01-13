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
    temp = []
    infile = open(filename, 'r')

    lines = [line.strip() for line in infile.readlines()]

    for line in lines:
        x.append(float(line.split(',')[0]))
        y.append(float(line.split(',')[1]))
        temp.append(float(line.split(',')[2]))

    return np.array(x), np.array(y), np.array(temp)

time = str(sys.argv[1])
filename = time + "_sec.dat"

x,y,temp = ParseFile(filename)
xi = np.linspace(min(x), max(x))
yi = np.linspace(min(y), max(y))
zi = griddata(x, y, temp, xi, yi, interp="linear")
x_grid, y_grid = np.meshgrid(xi,yi)

plotfile = time + "_sec_heat-eq_plot.pdf"
plotfile2 = time + "_sec_heat-eq_colorplot.pdf"
plotfile3 = time + "_sec_heat-eq_imshow.pdf"

#plotfile = time.replace('.','-') + "_sec_heat-eq_plot.pdf"
#plotfile2 = time.replace('.','-') + "_sec_heat-eq_colorplot.pdf"

matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'

plt.clf()
CS = plt.contour(xi, yi, zi)
plt.clabel(CS, inline=1, fontsize=10)
plt.title("Heated Bar (t = " + time + " seconds)")
plt.xlabel("x")
plt.ylabel("y")
#plt.show()
plt.savefig(plotfile, format="pdf")


plt.clf()
CS = plt.contour(xi,yi,zi,15,linewidths=0.5,colors='k')
CS = plt.contourf(xi,yi,zi,15,cmap=plt.cm.jet)
plt.colorbar()
plt.xlabel("x (m)")
plt.ylabel("y (m)")
plt.title("Heated Bar (t = " + time + " seconds)")
plt.clim(0,40)
#plt.show()
plt.savefig(plotfile2, format="pdf")

plt.clf()
zi = scipy.interpolate.griddata((x,y),temp,(x_grid,y_grid),method="linear")
plt.imshow(zi, vmin=temp.min(), vmax=temp.max(), origin="lower", extent=[x.min(), x.max(), y.min(), y.max()])
plt.xlabel("x (m)")
plt.ylabel("y (m)")
plt.title("Heated Bar (t = " + time + " seconds)")
plt.colorbar()
#plt.show()
plt.savefig(plotfile3, format="pdf")
