#!/usr/bin/python

from __future__ import print_function, division
import math
import numpy as np
import matplotlib.pyplot as plt
import sys

def ParseFile(filename):
    x = []
    u = []
    infile = open(filename, 'r')

    lines = [line.strip() for line in infile.readlines()]

    for line in lines:
        x.append(line.split(',')[0])
        u.append(line.split(',')[1])

    return x,u

time = str(sys.argv[1])
filename = time + "_sec.dat"

x,u = ParseFile(filename)

plotfile = time + "_sec_plot_CN2.pdf"

plt.clf()
plt.xlabel("Distance (m)")
plt.ylabel("Velocity (m/s)")
plt.ylim([-5,40])
plt.title("Infinite Parallel Plates (t = " + time + " seconds)")
plt.plot(x,u)
plt.savefig(plotfile, format="pdf")
