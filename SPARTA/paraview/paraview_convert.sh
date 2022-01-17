#!/bin/sh
cd surf
rm -rf *
pvpython ../../tools/surf2paraview.py ../../setup/data/data.MAV_stage_2 MAV_stage_2 
cd ../grid

rm -rf vals_MAV_stage_2_115km_0 
rm -rf vals_MAV_stage_2_115km_0.pvd 
echo 'Converting results of MAV_stage_2 at 115km (refinement 0) to ParaView...'
pvpython ../../tools/grid2paraview_original.py def/grid.MAV_stage_2_115km_0 vals_MAV_stage_2_115km_0 -r ../../setup/results_sparta/MAV_stage_2/vals_115km_0.*.dat 

rm -rf vals_MAV_stage_2_115km_1 
rm -rf vals_MAV_stage_2_115km_1.pvd 
echo 'Converting results of MAV_stage_2 at 115km (refinement 1) to ParaView...'
pvpython ../../tools/grid2paraview_original.py def/grid.MAV_stage_2_115km_1 vals_MAV_stage_2_115km_1 -r ../../setup/results_sparta/MAV_stage_2/vals_115km_1.*.dat 

rm -rf vals_MAV_stage_2_115km_2 
rm -rf vals_MAV_stage_2_115km_2.pvd 
echo 'Converting results of MAV_stage_2 at 115km (refinement 2) to ParaView...'
pvpython ../../tools/grid2paraview_original.py def/grid.MAV_stage_2_115km_2 vals_MAV_stage_2_115km_2 -r ../../setup/results_sparta/MAV_stage_2/vals_115km_2.*.dat 

rm -rf vals_MAV_stage_2_115km_3 
rm -rf vals_MAV_stage_2_115km_3.pvd 
echo 'Converting results of MAV_stage_2 at 115km (refinement 3) to ParaView...'
pvpython ../../tools/grid2paraview_original.py def/grid.MAV_stage_2_115km_3 vals_MAV_stage_2_115km_3 -r ../../setup/results_sparta/MAV_stage_2/vals_115km_3.*.dat 
