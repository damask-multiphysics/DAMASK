#!/usr/bin/env python

# renders the postprocessing routines acessible from everywhere.
# you need a ~/bin directory in your home folder
# if necessary, add the bin directory to your path by
# adding the following lines to your .bashrc file: 
#  PATH="$PATH:~/bin"
#  export PATH 

import os,sys,glob

bin_link = { \
				'pre' : [
							'marc_addUserOutput.py',
							'mentat_pbcOnBoxMesh.py',
							'mentat_spectralBox.py',
							'patchFromReconstructedBoundaries.py',
							'spectral_geomCheck.py',
							'spectral_minimalSurface.py',
							'spectral_vicinityOffset.py',
							'voronoi_randomSeeding.py',
							'voronoi_tessellation.py',
							],
				'post' : [
							'3Dvisualize.py',
							'addCauchy.py',
							'addDeterminant.py',
							'addDivergence.py',
							'addMises.py',
							'addNorm.py',
							'addStrainTensors.py',
              'addCompatibilityMismatch.py',
							'averageDown.py',
							'mentat_colorMap.py',
							'postResults.py',
							'spectral_iterationCount.py',
							'spectral_convergence.py',
							'spectral_scaleGeom.py',
							],
						}

compile = { \
				'pre' : [
							'voronoi_randomSeeding.f90',
							'voronoi_tessellation.f90',
							],
				'post' : [
							],
						}

execute = { \
				'pre' : [
							],
				'post' : [
							'make_postprocessingMath',
							],
						}

#homedir = os.getenv('HOME')
#basedir = os.path.join(os.path.dirname(sys.argv[0]),'..')
homedir = os.getenv('DAMASK_ROOT')
basedir = os.getenv('DAMASK_ROOT')+'/processing/'

for dir in bin_link:
  for file in bin_link[dir]:
    src = os.path.abspath(os.path.join(basedir,dir,file))
    if (file == ''):
      dst = os.path.abspath(os.path.join(homedir,'bin',dir))
    else:
      dst = os.path.abspath(os.path.join(homedir,'bin',os.path.splitext(file)[0]))
    print src,'-->',dst
    if os.path.lexists(dst):
      os.remove(dst)
    os.symlink(src,dst)
    
for dir in compile:
  for file in compile[dir]:
    src = os.path.abspath(os.path.join(basedir,dir,file))
    if os.path.isfile(src):
      print file
      os.system('rm %s.exe'%(os.path.splitext(src)[0]))
      os.system('ifort -O3 -parallel -o%s %s'%(os.path.splitext(src)[0],src))
 
modules = glob.glob('*.mod')
for module in modules:
  print module
  os.remove(module)

for dir in execute:
  for file in execute[dir]:
    cmd = os.path.abspath(os.path.join(basedir,dir,file))
    if os.access(cmd,os.X_OK):
      print cmd
      os.system('%s %s'%(cmd,os.path.join(basedir,dir)))
