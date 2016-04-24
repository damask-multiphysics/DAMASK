#!/bin/bash

for seeds in "$@"
do
  vtk_pointcloud $seeds

  vtk_addPointCloudData $seeds \
    --scalar microstructure,weight \
    --inplace \
    --vtk ${seeds%.*}.vtp \
    
done
