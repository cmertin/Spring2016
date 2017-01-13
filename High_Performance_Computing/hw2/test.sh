#!/bin/bash

million=1000000
tenMillion=10000000
hundredMillion=100000000
billion=1000000000
i=$million
qsort="quicksort"
scan="genericScan"

if [ ! -f $qsort ]; then
    echo "$qsort not found. Making..."
    make sort
fi

if [ ! -f $scan ]; then
    echo "$scan not found. Making..."
    make scan
fi

# Tests quicksort
echo -e "\n\t\tTESTING QUICKSORT"
while [ $i -le $billion ]; do
    title="Running for N = $i"
    size=${#title}
    echo ""
    echo -ne "\t"
    for ((j=0; j<$size; j++)); do echo -n '#'; done
    echo ""
    echo -e "\t$title"
    echo -ne "\t"
    for ((j=0; j<$size; j++)); do echo -n '#'; done
    echo ""
    ./quicksort $i
    i=$((i*10))
done

# Tests genericScan
i=$million
echo -e "\n\t\tTESTING GENERICSCAN"
while [ $i -le $billion ]; do
    title="Running for N = $i"
    size=${#title}
    echo ""
    echo -ne "\t"
    for ((j=0; j<$size; j++)); do echo -n '#'; done
    echo ""
    echo -e "\t$title"
    echo -ne "\t"
    for ((j=0; j<$size; j++)); do echo -n '#'; done
    echo ""
    ./genericScan $i
    i=$((i*10))
done
