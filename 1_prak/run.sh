#!/bin/sh


mpic++ -o MPI_RUN --std=c++11 MPI_calcul.cpp 

for i in $(seq 1 7)
do
mpirun -np $i ./MPI_RUN  $1 $2 
done