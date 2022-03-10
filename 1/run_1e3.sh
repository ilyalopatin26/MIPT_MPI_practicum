#!/bin/sh

mpic++ -o MPI_TIME  MPI_first_prak.cpp 

for i in $(seq 1 7)
do
 mpirun ./MPI_TIME  -np $i -N1000 
done