#!/bin/bash

#declare -a vel=("10" "20" "30" "40")
declare -a vel=("20" "30")

#for U in "${vel[@]}"; do
#    for sep in `seq 1 2 11`; do
#	./separation $U 10 38 15 15 $sep
#    done
#    echo "Finished U=$U from 1...11" | mail -s "U=$U Complete (First Half)" cmertin@cs.utah.edu #`seq 1 2 11`
    #echo "Finished U=$U from 13...23" | mail -s "U=$U Complete (Second Half)" cmertin@cs.utah.edu #`seq 13 2 23`
#done

# U startX startY width height separation
# baseGround is hardcoded to 18
for y in `seq 19 2 44`; do
    gnd=$((y-18))
    ./ground 10 20 $y 15 15 11
    #echo "Finished Ground Clearance of $gnd" | mail -s "Finished startY = $y" cmertin@cs.utah.edu
done

echo "Finished Plots" | mail -s "All Clearance Plots Finished" cmertin@cs.utah.edu
