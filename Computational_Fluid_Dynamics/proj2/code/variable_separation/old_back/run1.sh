#!/bin/bash

# Changing the separation by increments of 5
for i in `seq 1 5 26`; do
    sep=$i
    #u startX startY width height separation
    ./proj2 10 10 38 15 15 $sep
done
