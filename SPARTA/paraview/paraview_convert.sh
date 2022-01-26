#!/bin/sh
cd surf
rm -rf *
pvpython ../../tools/surf2paraview.py ../../setup/data/data.MAV_stage_2 MAV_stage_2 
cd ../grid

rm -rf vals_MAV_stage_2_100km_0 
rm -rf vals_MAV_stage_2_100km_0.pvd 
echo 'Converting results of MAV_stage_2 at 100km (refinement 0) to ParaView...'
pvpython ../../tools/grid2paraview_original.py def/grid.MAV_stage_2_100km_0 vals_MAV_stage_2_100km_0 -r ../../setup/results_sparta/MAV_stage_2/vals_100km_0.*.dat 

rm -rf vals_MAV_stage_2_100km_1 
rm -rf vals_MAV_stage_2_100km_1.pvd 
echo 'Converting results of MAV_stage_2 at 100km (refinement 1) to ParaView...'
pvpython ../../tools/grid2paraview_original.py def/grid.MAV_stage_2_100km_1 vals_MAV_stage_2_100km_1 -r ../../setup/results_sparta/MAV_stage_2/vals_100km_1.*.dat 

rm -rf vals_MAV_stage_2_150km_0 
rm -rf vals_MAV_stage_2_150km_0.pvd 
echo 'Converting results of MAV_stage_2 at 150km (refinement 0) to ParaView...'
pvpython ../../tools/grid2paraview_original.py def/grid.MAV_stage_2_150km_0 vals_MAV_stage_2_150km_0 -r ../../setup/results_sparta/MAV_stage_2/vals_150km_0.*.dat 

rm -rf vals_MAV_stage_2_150km_1 
rm -rf vals_MAV_stage_2_150km_1.pvd 
echo 'Converting results of MAV_stage_2 at 150km (refinement 1) to ParaView...'
pvpython ../../tools/grid2paraview_original.py def/grid.MAV_stage_2_150km_1 vals_MAV_stage_2_150km_1 -r ../../setup/results_sparta/MAV_stage_2/vals_150km_1.*.dat 

rm -rf vals_MAV_stage_2_200km_0 
rm -rf vals_MAV_stage_2_200km_0.pvd 
echo 'Converting results of MAV_stage_2 at 200km (refinement 0) to ParaView...'
pvpython ../../tools/grid2paraview_original.py def/grid.MAV_stage_2_200km_0 vals_MAV_stage_2_200km_0 -r ../../setup/results_sparta/MAV_stage_2/vals_200km_0.*.dat 

rm -rf vals_MAV_stage_2_200km_1 
rm -rf vals_MAV_stage_2_200km_1.pvd 
echo 'Converting results of MAV_stage_2 at 200km (refinement 1) to ParaView...'
pvpython ../../tools/grid2paraview_original.py def/grid.MAV_stage_2_200km_1 vals_MAV_stage_2_200km_1 -r ../../setup/results_sparta/MAV_stage_2/vals_200km_1.*.dat 

rm -rf vals_MAV_stage_2_250km_0 
rm -rf vals_MAV_stage_2_250km_0.pvd 
echo 'Converting results of MAV_stage_2 at 250km (refinement 0) to ParaView...'
pvpython ../../tools/grid2paraview_original.py def/grid.MAV_stage_2_250km_0 vals_MAV_stage_2_250km_0 -r ../../setup/results_sparta/MAV_stage_2/vals_250km_0.*.dat 

rm -rf vals_MAV_stage_2_250km_1 
rm -rf vals_MAV_stage_2_250km_1.pvd 
echo 'Converting results of MAV_stage_2 at 250km (refinement 1) to ParaView...'
pvpython ../../tools/grid2paraview_original.py def/grid.MAV_stage_2_250km_1 vals_MAV_stage_2_250km_1 -r ../../setup/results_sparta/MAV_stage_2/vals_250km_1.*.dat 

rm -rf vals_MAV_stage_2_300km_0 
rm -rf vals_MAV_stage_2_300km_0.pvd 
echo 'Converting results of MAV_stage_2 at 300km (refinement 0) to ParaView...'
pvpython ../../tools/grid2paraview_original.py def/grid.MAV_stage_2_300km_0 vals_MAV_stage_2_300km_0 -r ../../setup/results_sparta/MAV_stage_2/vals_300km_0.*.dat 

rm -rf vals_MAV_stage_2_300km_1 
rm -rf vals_MAV_stage_2_300km_1.pvd 
echo 'Converting results of MAV_stage_2 at 300km (refinement 1) to ParaView...'
pvpython ../../tools/grid2paraview_original.py def/grid.MAV_stage_2_300km_1 vals_MAV_stage_2_300km_1 -r ../../setup/results_sparta/MAV_stage_2/vals_300km_1.*.dat 
