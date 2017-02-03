#!/usr/bin/env bash

if [[ "$1" == "-f" || "$1" == "--float" ]]; then
  shift
  arg='--float'
else
  arg=''
fi 

for geom in "$@"
do
  geom_toTable $arg \
  < $geom \
  | \
  vtk_rectilinearGrid > ${geom%.*}.vtk

  geom_toTable $arg \
  < $geom \
  | \
  vtk_addRectilinearGridData \
    --data microstructure \
    --inplace \
    --vtk ${geom%.*}.vtk
  rm ${geom%.*}.vtk
done
