#!/bin/sh
mpirun -np 16 spa_ < in.MAV_stage_2_100km | tee ../results_sparta/MAV_stage_2/stats_100km.dat
mpirun -np 16 spa_ < in.MAV_stage_2_150km | tee ../results_sparta/MAV_stage_2/stats_150km.dat
mpirun -np 16 spa_ < in.MAV_stage_2_200km | tee ../results_sparta/MAV_stage_2/stats_200km.dat
mpirun -np 16 spa_ < in.MAV_stage_2_250km | tee ../results_sparta/MAV_stage_2/stats_250km.dat
mpirun -np 16 spa_ < in.MAV_stage_2_300km | tee ../results_sparta/MAV_stage_2/stats_300km.dat
