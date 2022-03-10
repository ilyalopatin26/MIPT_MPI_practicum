#!/bin/sh

mpic++ -o MPI_TIME --std=c++11  timing.cpp 

for i in $(seq 1 7)
do
 mpirun ./MPI_TIME  -np $i -N$1
done