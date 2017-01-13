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
        if "Separation" in line:
            separation = int(line.split(' ')[2])
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
        if "Width" in line:
            width = int(line.split(' ')[2])
        if "Nx" in line:
            Nx = int(line.split(' ')[2])
        if not line.startswith("#"):
            x.append(float(line.split(',')[0]))
            y.append(float(line.split(',')[1]))
            u.append(float(line.split(',')[2]))
            v.append(float(line.split(',')[3]))
            P.append(float(line.split(',')[4]))
            omega.append(float(line.split(',')[5]))

    return startX,startY,separation,dx,dy,height,width,Nx,np.array(x),np.array(y),np.array(u),np.array(v),np.array(P),np.array(omega)

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

    return ptotal*dx*dy
        

#filename = "train_10_15_15_0_data.dat"
filename = "train_10_15_15_9_data.dat"

startX,startY,separation,dx,dy,height,width,Nx,x,y,u,v,P,omega = ParseFile(filename)

v = np.reshape(v, (-1, Nx))
u = np.reshape(u, (-1, Nx))
P = np.reshape(P, (-1, Nx))

midPoint = int(math.ceil(height/2)) + startY

cv_width = range(startX-20,startX-2)
cv_width.reverse()

#print(startX + height)

#print(v[startX, startY+height])

width_do = 17

alpha = CalcAlpha(u, v, dx, dy, startX, startY, height, width_do)

#print("alpha: " + str(alpha))

P = CalcPressure(P, startX, startY, dx, dy, height, width_do)


U = 10
base_filename = "train_" + str(U) + "_15_15_"
velocities = range(5,101,5)
velocities_true = range(0,101)
sep_list = range(1,56,2)
fd_true = []
error = 0.033
fd = []
fd_err = []
fd_sub = []
fd_sub_err = []
fd_ratio = []
fd_ratio_err = []
sep_final = []
sep_ratio = []
print(sep_list)

for sep in sep_list:
    if sep > 25:
        break
    
    filename = base_filename + str(sep) + "_data.dat"
    print("Sep " + str(sep) + ": ", end="")
    startX,startY,separation,dx,dy,height,width,Nx,x,y,u,v,P,omega = ParseFile(filename)
    exact = 0.5 * 2.2 * U * U * height * dy
    startX = startX + width + separation
    v = np.reshape(v, (-1, Nx))
    u = np.reshape(u, (-1, Nx))
    P = np.reshape(P, (-1, Nx))
    cv_width = range(startX-20, startX-2)
    cv_width.reverse()
    width_do = 17

    sep_final.append(sep)
    sep_ratio.append(sep/width)
    alpha = CalcAlpha(u, v, dx, dy, startX, startY, height, width_do)
    P = CalcPressure(P, startX, startY, dx, dy, height, width_do)
    fd.append(P + alpha)
    fd_err.append((P + alpha)*error)
    fd_sub.append(P + alpha - exact)
    fd_sub_err.append((P + alpha - exact)*error)
    fd_ratio.append((P + alpha)/exact)
    fd_ratio_err.append((P + alpha)/exact * error)
    error_d = "%.5f" % (error * (P + alpha))
    error_dd = "%.5f" % (error * (P + alpha - exact))
    drag = "%.4f" % (P + alpha)
    d_drag = "%.4f" % (P + alpha - exact)
    print(drag + " +/- " + error_d +  " (", end="")
    print(d_drag + " +/- " + error_dd + ")")

print("Printing Style for PGF Plots")
print("============================\n\n")
print("\\begin{tikzpicture}")
print("\\begin{axis}[xlabel={$\\frac{Separation}{Width}$}, ylabel={$\Delta F_{d}$}, title = {Change in $F_{d}$}]")#, xtick={",end="")
#for i in xrange(0, len(sep_ratio)):
#    x = "%.4f" % sep_ratio[i]
#    print(x + ",",end="")
#print("}, xticklabels={",end="")
#for i in xrange(0, len(sep_final)):
#    print("$\\frac{" + str(sep_final[i]) + "}{" + str(width) + "}$",end="")
#    if(i != len(sep_final)-1):
#        print(",",end="")
#print("}]")
print("\\addplot+[error bars/.cd, y dir=both, y explicit] coordinates{")#color=blue, mark=*, thick] coordinates{")#, end="")
for i in xrange(0, len(fd_sub)):
    x = "%.4f" % sep_ratio[i]
    y = "%.4f" % fd_sub[i]
    err = "%.4f" % fd_sub_err[i]
    print("(" + x + "," + y + ") +- (" + err + "," + err + ")")
print("};")
print("\\addlegendentry{$U = " +  str(U) + "$}")
print("\\end{axis}")
print("\\end{tikzpicture}")

print("\n\nRatio")
print("=====\n\n")
print("\\begin{tikzpicture}")
print("\\begin{axis}[xlabel={$\\frac{Separation}{Width}$}, ylabel={$\\frac{F_{d}^{(car)}}{F_{d}}$}, title={Ratio of Train Car to Single Block}, y tick label style={/pgf/number format/.cd,fixed,fixed zerofill,precision=3,/tikz/.cd}]")#, xtick={",end="")
#for i in xrange(0, len(sep_ratio)):
#    x = "%.4f" % sep_ratio[i]
#    print(x + ",",end="")
#print("}, xticklabels={",end="")
#for i in xrange(0, len(sep_final)):
#    print("$\\frac{" + str(sep_final[i]) + "}{" + str(width) + "}$",end="")
#    if(i != len(sep_final)-1):
#        print(",",end="")
#print("}]")
print("\\addplot[mark=*, thick, color=blue] coordinates{")#color=blue, mark=*, thick] coordinates{")#, end="")
for i in xrange(0, len(fd_ratio)):
    x = "%.4f" % sep_ratio[i]
    y = "%.4f" % fd_ratio[i]
    err = "%.4f" % fd_ratio_err[i]
    print("(" + x + "," + y + ")")
print("};")
print("\\addlegendentry{$U = " +  str(U) + "$}")
print("\\end{axis}")
print("\\end{tikzpicture}")
'''
for vel in velocities_true:
    fd_temp = 0.5 * 2.2 * vel * vel * height * dy
    fd_true.append(fd_temp)

for vel in velocities:
    filename = "train_" + str(vel) + base_filename
    startX,startY,dx,dy,height,Nx,x,y,u,v,P,omega = ParseFile(filename)

    v = np.reshape(v, (-1, Nx))
    u = np.reshape(u, (-1, Nx))
    P = np.reshape(P, (-1, Nx))

    midPoint = int(math.ceil(height/2)) + startY

    cv_width = range(startX-20,startX-2)
    cv_width.reverse()
    width_do = 17

    alpha = CalcAlpha(u, v, dx, dy, startX, startY, height, width_do)

    P = CalcPressure(P, startX, startY, dx, dy, height, width_do)

    fd.append(P + alpha)
    fd_true_temp = 0.5 * 2.2 * vel * vel * height * dy
    temp_error = abs(P + alpha - fd_true_temp)/fd_true_temp * 100
    error.append(temp_error)



plt.clf()
plt.rc("text", usetex=True)
plt.rc("font", family="serif")
plt.plot(velocities_true, fd_true, 'k', label="$Exact$")
plt.plot(velocities, fd, "bo", label="$Control\ Volume$")
plt.xlabel("$Velocity\ (u)$")
plt.ylabel("$Pressure\ (P)$")
plt.title("Pressure vs Velocity Using Exact Equation and Control Volumes")
plt.legend(loc=2)
plt.savefig("Pressure_vs_Vel.pdf", filetype="pdf", bbox_inches="tight")#, pad_inches=0)
#plt.show()

plt.clf()
plt.rc("text", usetex=True)
plt.rc("font", family="serif")
plt.plot(velocities, error, 'k', color="red")
plt.xlabel("$Velocity\ (u)$")
plt.ylabel("$\%\ Error$")
plt.title("Percent Error For Velocities")
#plt.show()
plt.savefig("Percent_Error.pdf", filetype="pdf", bbox_inches="tight")#, pad_inches=0)

print("Printing style for PGF Plots")
print("============================\n\n")
print("\\begin{tikzpicture}")
print("\\begin{axis}[xlabel={$Velocity\ (U)$}, ylabel={$Pressure\ (P)$}, legend pos=north west, title={Pressure vs Velocity for Exact and Control Volume}]")
print("\\addplot[thick, color=black] coordinates{",end="")
for i in xrange(0, len(velocities_true)):
    fd_str = "%.3f" % fd_true[i]
    print("(" + str(velocities_true[i]) + ',' + fd_str + ")",end="")
print("};")
print("\\addlegendentry{Exact}")

print("\\addplot[color=blue, mark=*, only marks] coordinates{",end="")
for i in xrange(0, len(velocities)):
    fd_str = "%.3f" % fd[i]
    print("(" + str(velocities[i]) + "," + fd_str + ")",end="")
print("};")
print("\\addlegendentry{Control Volume}")
print("\\end{axis}")
print("\\end{tikzpicture}")


print("\n\n")
print("\\begin{tikzpicture}")
print("\\begin{axis}[xlabel={$Velocity\ (U)$}, ylabel={$\%\ Error$}, legend pos=north east, title={Percent Error for $F_{d}$}, y tick label style={/pgf/number format/.cd,fixed,fixed zerofill,precision=3,/tikz/.cd}]")
print("\\addplot[thick, color=red] coordinates{", end="")
for i in xrange(0, len(velocities)):
    err_str = "%.4f" % error[i]
    print("(" + str(velocities[i]) + "," + err_str + ")", end="")
print("};")
print("\\end{axis}")
print("\\end{tikzpicture}")
'''