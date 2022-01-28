#!/bin/sh
mpirun -np 16 spa_ < in.MAV_stage_2_100km | tee ../results_sparta/MAV_stage_2/stats_100km.dat
mpirun -np 16 spa_ < in.MAV_stage_2_125km | tee ../results_sparta/MAV_stage_2/stats_125km.dat
mpirun -np 16 spa_ < in.MAV_stage_2_150km | tee ../results_sparta/MAV_stage_2/stats_150km.dat
mpirun -np 16 spa_ < in.MAV_stage_2_175km | tee ../results_sparta/MAV_stage_2/stats_175km.dat
mpirun -np 16 spa_ < in.MAV_stage_2_200km | tee ../results_sparta/MAV_stage_2/stats_200km.dat
mpirun -np 16 spa_ < in.MAV_stage_2_225km | tee ../results_sparta/MAV_stage_2/stats_225km.dat
mpirun -np 16 spa_ < in.MAV_stage_2_250km | tee ../results_sparta/MAV_stage_2/stats_250km.dat
mpirun -np 16 spa_ < in.MAV_stage_2_275km | tee ../results_sparta/MAV_stage_2/stats_275km.dat
mpirun -np 16 spa_ < in.MAV_stage_2_300km | tee ../results_sparta/MAV_stage_2/stats_300km.dat
mpirun -np 16 spa_ < in.MAV_stage_2_350km | tee ../results_sparta/MAV_stage_2/stats_350km.dat
mpirun -np 16 spa_ < in.MAV_stage_2_400km | tee ../results_sparta/MAV_stage_2/stats_400km.dat
mpirun -np 16 spa_ < in.MAV_stage_2_450km | tee ../results_sparta/MAV_stage_2/stats_450km.dat
mpirun -np 16 spa_ < in.MAV_stage_2_500km | tee ../results_sparta/MAV_stage_2/stats_500km.dat
