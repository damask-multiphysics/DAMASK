#!/bin/bash

for geom in "$@"
do
  vtk_rectilinearGrid \
    --geom $geom

  geom_toTable \
  < $geom \
  | \
  vtk_addRectilinearGridData \
    --scalar microstructure \
    --inplace \
    --vtk ${geom%.*}.vtr
done
