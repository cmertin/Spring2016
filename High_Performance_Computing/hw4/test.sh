#!/bin/bash

for i in `seq 1 10`; do
    mpirun -np 16 ./sort 1
done
