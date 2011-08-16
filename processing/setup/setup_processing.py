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
							'marc_addUserOutput',
							'mentat_pbcOnBoxMesh',
							'mentat_spectralBox',
							'patchFromReconstructedBoundaries',
							'spectral_geomCheck',
							'voronoi_randomSeeding',
							'voronoi_tessellation',
							],
				'post' : [
							'3Dvisualize',
							'addCauchy',
							'addDeterminant',
							'addDivergence',
							'addMises',
							'addNorm',
							'addStrainTensors',
							'averageDown',
							'mentat_colorMap',
							'postResults',
							'spectral_iterationCount',
							'spectral_convergence',
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

homedir = os.getenv('HOME')
basedir = os.path.join(os.path.dirname(sys.argv[0]),'..')
for dir in bin_link:
  for file in bin_link[dir]:
    src = os.path.abspath(os.path.join(basedir,dir,file))
    if (file == ''):
      dst = os.path.abspath(os.path.join(homedir,'bin',dir))
    else:
      dst = os.path.abspath(os.path.join(homedir,'bin',file))
    print src,'-->',dst
    if os.path.lexists(dst):
      os.remove(dst)
    os.symlink(src,dst)

for dir in compile:
  for file in compile[dir]:
    src = os.path.abspath(os.path.join(basedir,dir,file))
    if os.path.isfile(src):
      print file
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
