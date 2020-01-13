#!/usr/bin/env bash

for seeds in "$@"
do
  vtk_pointCloud $seeds

  vtk_addPointCloudData $seeds \
    --data microstructure,weight \
    --vtk ${seeds%.*}.vtp \
    
done
