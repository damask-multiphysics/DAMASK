#!/usr/bin/env bash

for geom in $(ls geom/*.geom)
do
  base=${geom%.geom}
  base=${base#geom/}
  name=${base}_thermal
  vtr=${base}.vtr

  [[ -f ${name}.spectralOut ]] || \
  DAMASK_spectral \
    --workingdir ./ \
    --load thermal.load \
    --geom $geom \
    > ${name}.out
  
  if [ ! -f postProc/${name}_inc10.txt ]
  then
    postResults ${name}.spectralOut \
      --ho temperature \
      --cr f,fe,fi,fp,p \
      --split \
      --separation x,y,z \

    addCauchy postProc/${name}_inc*.txt \

    addDeviator postProc/${name}_inc*.txt \
      --spherical \
      --tensor p,Cauchy \

    addDisplacement postProc/${name}_inc*.txt \
      --nodal \

  fi

  geom_check ${geom}
  
  for inc in {00..10}
  do
    echo "generating postProc/${name}_inc${inc}.vtr"
     cp geom/${vtr} postProc/${name}_inc${inc}.vtr
     vtk_addRectilinearGridData \
       postProc/${name}_inc${inc}.txt \
       --vtk postProc/${name}_inc${inc}.vtr \
       --data 'sph(p)','sph(Cauchy)',temperature \
       --tensor f,fe,fi,fp,p,Cauchy \
      
    vtk_addRectilinearGridData \
      postProc/${name}_inc${inc}_nodal.txt \
      --vtk postProc/${name}_inc${inc}.vtr \
      --data 'avg(f).pos','fluct(f).pos' \

  done
done
