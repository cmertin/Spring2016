#!/bin/bash

#prog    u startX startY width height separation

j=1

for i in `seq 5 5 100`; do
    ./ground $i 38 38 15 15 0
    #echo "Finished U=$i for the case of the single block" | mail -s "Run $j/19 Complete" cmertin@cs.utah.edu
    echo "Finished U=$i"
    j=$((j+1))	
done

#echo "test" | mail -s "test" cmertin@cs.utah.edu


#
