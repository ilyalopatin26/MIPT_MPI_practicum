#!/bin/sh

mpic++ -o MPI_FP  timing.cpp 

for i in $(seq 1 7)
do
 mpirun ./MPI_FP  -np $i -N1000 
done