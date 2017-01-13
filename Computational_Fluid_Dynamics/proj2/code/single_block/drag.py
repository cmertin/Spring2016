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
    dx = None
    dy = None
    startX = None
    startY = None

    infile = open(filename, 'r')

    lines = [line.strip() for line in infile.readlines()]

    for line in lines:
        if "Start X" in line:
            startX = int(line.split(' ')[3])
        if "Start Y" in line:
            startY = int(line.split(' ')[3])
        if "dy" in line:
            dy = float(line.split(' ')[2])
        if "dx" in line:
            dx = float(line.split(' ')[2])
        if "Height" in line:
            height = int(line.split(' ')[2])
        if "Nx" in line:
            Nx = int(line.split(' ')[2])
        if not line.startswith("#"):
            x.append(float(line.split(',')[0]))
            y.append(float(line.split(',')[1]))
            u.append(float(line.split(',')[2]))
            v.append(float(line.split(',')[3]))
            P.append(float(line.split(',')[4]))
            omega.append(float(line.split(',')[5]))

    return startX,startY,dx,dy,height,Nx,np.array(x),np.array(y),np.array(u),np.array(v),np.array(P),np.array(omega)

def CalcAlpha(u, v, dx, dy, startX, startY, height, cv_width):
    top = 0
    bottom = 0
    left = 0
    right = 0

    startY = startY - cv_width
    endY = startY + height + 2 * cv_width
    startX = startX - cv_width
    endX = startX + height + 2 * cv_width
    for i in xrange(startY, endY):
        ur = u[startX, i]
        vr = v[startX, i]
        ul = u[startX-cv_width, i]
        vl = v[startX-cv_width, i]
        right += ur*dx*(ur + vr)
        left += ul*dx*(ul+vl)

    for i in xrange(startX, endX):
        vt = v[i, startY+height]
        vb = v[i, startY]
        ut = u[i, startY+height]
        ub = u[i, startY]
        top += vt*dy*(vt + ut)
        bottom += vb*dy*(vb + ub)

    alpha = right - left + top - bottom
    return alpha

def CalcPressure(P, startX, startY, dx, dy, height, cv_width):
    ptotal = 0
    startY = startY - cv_width
    endY = startY + height + 2 * cv_width
    startX = startX - cv_width
    endX = startX + height + 2 * cv_width
    for j in xrange(startY, endY):
        for i in xrange(startX, endX):
            ptotal += P[i][j]
    print(ptotal*dx*dy)
        

#filename = "train_10_15_15_0_data.dat"
filename = "train_20_15_15_0_data.dat"

startX,startY,dx,dy,height,Nx,x,y,u,v,P,omega = ParseFile(filename)

v = np.reshape(v, (-1, Nx))
u = np.reshape(u, (-1, Nx))
P = np.reshape(P, (-1, Nx))

midPoint = int(math.ceil(height/2)) + startY

cv_width = range(startX-20,startX-2)
cv_width.reverse()

print(startX + height)

print(v[startX, startY+height])

width_do = 17

alpha = CalcAlpha(u, v, dx, dy, startX, startY, height, width_do)

print("alpha: " + str(alpha))

CalcPressure(P, startX, startY, dx, dy, height, width_do)

