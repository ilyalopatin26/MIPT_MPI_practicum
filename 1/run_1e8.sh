#!/bin/sh

mpic++ -o MPI_FP  MPI_first_prak.cpp 

for i in $(seq 1 7)
do
 mpirun ./MPI_FP  -np $i -N100000000 
done