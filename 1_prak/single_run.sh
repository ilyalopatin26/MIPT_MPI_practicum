#!/bin/sh

mpic++ -o MPI_RUN --std=c++11 MPI_calcul.cpp 

mpirun -np $1 ./MPI_RUN  $2 $3
