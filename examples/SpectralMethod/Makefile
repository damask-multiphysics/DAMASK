
include ${PETSC_DIR}/lib/petsc-conf/variables
include ${PETSC_DIR}/lib/petsc-conf/rules


run16x16x16:
	-@${MPIEXEC} -n 2 DAMASK_spectral -l tensionX.load -g 20grains16x16x16.geom

run32x32x32:
	-@${MPIEXEC} -n 4 DAMASK_spectral -l tensionX.load -g 20grains32x32x32.geom

run64x64x64:
	-@${MPIEXEC} -n 8 DAMASK_spectral -l tensionX.load -g 20grains64x64x64.geom

