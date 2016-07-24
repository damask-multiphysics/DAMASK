#!/usr/bin/env bash

for seeds in "$@"
do
  vtk_pointcloud $seeds

  vtk_addPointcloudData $seeds \
    --scalar microstructure,weight \
    --inplace \
    --vtk ${seeds%.*}.vtp \
    
done
