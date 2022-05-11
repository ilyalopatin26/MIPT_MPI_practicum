#!/bin/sh


mpic++ -o MPI_RUN --std=c++11 PDE_MPI.cpp 

for i in $(seq 1 12)
do
mpirun -np $i ./MPI_RUN  -t 
done