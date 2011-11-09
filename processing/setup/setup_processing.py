#!/usr/bin/env python

# Makes postprocessing routines acessible from everywhere.

import os,sys,glob

import damask_tools; reload(damask_tools)
damask_tools.DAMASK_TOOLS().check_env()

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

damask_root = os.getenv('DAMASK_ROOT')
if damask_root == '' or damask_root == None: damask_root = os.path.join(os.path.dirname(sys.argv[0]),'../../')
damask_bin  = os.getenv('DAMASK_BIN')
if damask_bin  == '' or damask_bin  == None: damask_bin = os.path.join(damask_root,'bin/')
basedir = os.path.join(damask_root,'processing/')

for dir in bin_link:
  for file in bin_link[dir]:
    src = os.path.abspath(os.path.join(basedir,dir,file))
    if (file == ''):
      sym_link = os.path.abspath(os.path.join(damask_bin,dir))
    else:
      sym_link = os.path.abspath(os.path.join(damask_bin,os.path.splitext(file)[0]))
    print sym_link,'-->',src
    if os.path.lexists(sym_link):
      os.remove(sym_link)    
    os.symlink(src,sym_link)
    
    #--- uncomment next lines to remove your old symbolic links in ~/bin
    #old_link=sym_link.replace(damask_root,os.getenv('HOME'))
    #if os.path.lexists(old_link):
    #  os.remove(old_link)

    
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
