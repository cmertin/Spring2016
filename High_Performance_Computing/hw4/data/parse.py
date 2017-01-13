#!/usr/bin/python

from __future__ import print_function, division
import math
import numpy as np
import matplotlib.pyplot as plt
import sys

def Average(times, num):
    split = int(len(times)/num)
    data = [times[x:x+split] for x in range(0,len(times),split)]
    avg_data = []
    for item in data:
        avg_data.append(np.average(item))

    return avg_data

def ParseFile(filename, num):
    time = []
    infile = open(filename, 'r')

    lines = [line.strip() for line in infile.readlines()]

    for line in lines:
        if "Max Execution" in line:
            time.append(float(line.split()[3]))
    return Average(time, num)

# Outputs the Tikz code for p vs time
def Output_TeX2(nodes, x, y, shapes, colors):
    print("")
    print("Printing alternative tikz-pgfplot (time vs p)")
    print("")
    for i in xrange(0, len(x)):
        print("% n = " + str(x[i]))
        print("\\addplot[color=" + str(colors[i % len(colors)]) + ", mark=" + str(shapes[i % len(shapes)]) + ",]")
        print("\tcoordinates{",end="")
        for j in xrange(0, len(nodes)):
            print("(" + str(nodes[j]) + ", " + str(y[j][i]) + ")", end="")
        print("};")

# Outputs the Tikz code for n/p vs time
def Output_TeX(nodes, x, y, shapes, colors):
    print("")
    print("Printing for tikz-pgfplot (time vs n/p)")
    print("")
    for i in xrange(0, len(nodes)):
        print("% Nodes = " + str(nodes[i]))
        print("\\addplot[color=" + str(colors[i % len(colors)]) + ", mark=" + str(shapes[i % len(shapes)]) + ",]")
        print("\tcoordinates{",end="")
        for j in xrange(0, len(x)):
            nicex = str(x[j])
            nicey = "%.8f" % y[i][j]
            print('(' + nicex + "," + nicey + ')',end="")
        print("};")

# Outputs the body of a LaTeX table
def Output_TeXTable(nodes, data, vals):
    print("")
    print("Printing for LaTeX Table")
    print("")
    print("\\begin{tabular}{",end="")
    for i in xrange(0, len(vals)+1):
        print(" c",end="")
    print("}")
    print("\hline\hline")
    print("Nodes ", end="")
    for val in vals:
        print("& $10^{" + str(int(math.ceil(math.log(float(val),10)))) + "}$ ",end="")
    print("\\\\")
    print("\hline")
    for i in xrange(0,len(nodes)):
        print(str(nodes[i]),end="")
        for time in data[i]:
            nicetime = "%.4f" % time
            print(" & " + str(nicetime),end="")
        print("\\\\")
    print("\hline\hline")
    print("\end{tabular}")

def Output_Table(nodes, data, vals):
    print("")
    print("Printing for Table for Site")
    print("")
    labels = ["Nodes"]
    for val in vals:
        labels.append("10<sup>" + str(int(math.ceil(math.log(float(val),10)))) + "</sup>")
    for label in labels:
        print("|" + label, end="")
    print("|")
    for i in xrange(0,len(labels)):
        print("|:" + '='*(len(labels[i])-2) + ":",end="")
    print("|")
    for i in xrange(0,len(nodes)):
        print('|' + str(nodes[i]),end="")
        for time in data[i]:
            nicetime = "%.4f" % time
            print('|' + str(nicetime),end="")
        print("|")
        
vals = [1,10,100,1000,10000,100000,1000000,10000000]
nodes = [1,2,4,8,16,32]
tempvals = [1,10,100]
all_times = []
shapes = ["square*", "diamond*", "*", "triangle*"]
colors = ["blue", "orange", "red", "green", "gray", "violet", "teal", "olive", "lime", "magenta"]

for node in nodes:
    filename = "bitonic-" + str(node) + ".dat"
    times = ParseFile(filename, len(vals))
    all_times.append(times)

Output_Table(nodes, all_times, vals)
    
Output_TeXTable(nodes, all_times, vals)
    
Output_TeX(nodes, vals, all_times, shapes, colors)

Output_TeX2(nodes, vals, all_times, shapes, colors)

