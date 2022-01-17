#!/bin/sh
module load gnu/7
module load openmpi
mpirun -np 8 spa_ < in.MAV_stage_2_115km | tee ../results_sparta/MAV_stage_2/stats_115km.dat
