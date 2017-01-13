#!/bin/bash

for i in `seq 31 5 56`; do
    sep=$i
    #u startX startY width height separation
    ./proj2 10 10 38 15 15 $sep
done
