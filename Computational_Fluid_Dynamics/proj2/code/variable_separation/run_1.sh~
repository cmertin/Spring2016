#!/bin/bash

declare -a vel=("10" "20" "30" "40")

for U in "${vel[@]}"; do
    for sep in `seq 1 2 30`; do
	echo "./separation $U 10 38 15 15 $sep"
    done
    echo "Finished U=$U from 1...11" | mail -s "U=$U Complete (First Half)" cmertin@cs.utah.edu #`seq 1 2 11`
    #echo "Finished U=$U from 13...23" | mail -s "U=$U Complete (Second Half)" cmertin@cs.utah.edu #`seq 13 2 23`
done
