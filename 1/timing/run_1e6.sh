#!/bin/sh

mpic++ -o MPI_TIME  timing.cpp 

for i in $(seq 1 7)
do
 mpirun ./MPI_TIME  -np $i -N1000000 
done