#!/bin/bash

executable="train"

if [ ! -f $executable ]; then
    echo "Executable \"$executable\" not found. Attempted to make with Makefile"
    if [ ! -f Makefile ]; then
	echo "Makefile not found. Exiting"
	exit
    fi
    make
    echo "Make Complete. Continuing"
fi

declare -a vel=("10" "20" "30" "40")

for U in "${vel[@]}"; do
    for sep in `seq 1 2 25`; do
	./train $U 10 38 15 15 $sep
    done
    echo "Finished U=$U from 1...11" | mail -s "U=$U Complete" [EMAIL]
done
