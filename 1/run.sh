#!/bin/sh

mpic++ -o MPI_FP --std=c++11  MPI_first_prak.cpp 

for i in $(seq 1 7)
do
 mpirun ./MPI_FP  -np $i -N$1
done