#!/bin/sh
module load openmpi
export PATH="/mnt/c/Users/jerem/Downloads/sparta-master/build/src":$PATH
mpirun -np 8 spa_ < in.MAV_stage_2_115km | tee ../results_sparta/MAV_stage_2/stats_115km.dat
