#!/bin/bash

# Changing the separation by increments of 5
for i in `seq 1 2 30`; do
    sep=$i
    echo $i
    #u startX startY width height separation
    ./separation 40 10 38 15 15 $sep
    echo "Finished Separation $sep for the case of the single block" | mail -s "Run Complete" cmertin@cs.utah.edu
done