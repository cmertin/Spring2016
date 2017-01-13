#!/usr/bin/python

from __future__ import division, print_function
from math import sqrt

N = 50
nr = int((sqrt(N)-1)/2)
old = []
new = []

for i in xrange(1, nr+1):
    for j in xrange(1, nr+1):
        old.append((2 * (i - 1) + 1) * (2 * nr + 1) + 2 * j + 1)

for i in xrange(0, nr):
    for j in xrange(0, nr):
        new.append((2 * i + 1) * (2 * nr + 1) + 2 * (j + 1) + 1)

for i in xrange(0, len(old)):
    print(old[i], new[i])
