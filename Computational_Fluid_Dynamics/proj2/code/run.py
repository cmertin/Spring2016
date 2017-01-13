#!/usr/bin/python

from __future__ import print_function, division
import os
import os.path


# Check to see if the executable exists
executable = "./train"
if os.path.isfile(executable) == False:
    print("Executable not found. Attempting to use Makefile")
    if os.path.isfile("Makefile") == False:
        print("Makefile not found. Exiting.")
        exit()
    else:
        os.system("make")
    print("Make successful. Continuting...")

print("\n\t##############################")
print("\tNeed Parameters Before Running")
print("\t##############################\n")
print("\t===================")
print("\tDefault Parameters:")
print("\t===================")
print("\tRe = 100")
print("\t[min x, max x] = [0, 2.0]")
print("\t[min y, max y] = [0, 2.0]")
print("\tNx = 100")
print("\tNy = 100")
print("\n\t========================")
print("\tUser Defined Parameters:")
print("\t========================")
vel = str(int(raw_input("\tVelocity (U): ")))
stX = str(int(raw_input("\tIndex of Train to start at in X-direction [0,100]: ")))
stY = str(int(raw_input("\tIndex of Train to start at in Y-direction [0,100]: ")))
wdth = str(int(raw_input("\tUnits of train width: ")))
hght = str(int(raw_input("\tUnits of train height: ")))
sep = str(int(raw_input("\tUnits of separation between train cars: ")))

sys_cmd = executable + " " + vel + " " + stX + " " + stY + " " + wdth + " " + hght + " " + sep

os.system(sys_cmd)
