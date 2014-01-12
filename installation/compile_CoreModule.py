#!/usr/bin/env python

# Makes postprocessing routines acessible from everywhere.

import os,sys,glob,string,subprocess
from damask import Environment

damaskEnv = Environment()
baseDir = damaskEnv.relPath('installation/')
codeDir = damaskEnv.relPath('code/')

try:
  imklroot = damaskEnv.options['IMKLROOT']
except:
  imklroot = os.getenv('IMKLROOT')
  if imklroot == None: imklroot == ''           # env not set

try:
  acmlroot = damaskEnv.options['ACMLROOT']
except:
  acmlroot = os.getenv('ACMLROOT')
  if acmlroot == None: acmlroot == ''           # env not set

try:
  lapackroot = damaskEnv.options['LAPACKROOT']
except:
  lapackroot = os.getenv('LAPACKROOT')
  if lapackroot == None: lapackroot == ''       # env not set

try:
  fftwroot = damaskEnv.options['FFTWROOT']
except:
  fftwroot = os.getenv('FFTWROOT')
  if fftwroot == None: fftwroot == ''       # env not set

try:
  compiler = damaskEnv.options['F90']
except:
  compiler = os.getenv('F90')
  if compiler == None: compiler == ''       # env not set

compilers = ['ifort','gfortran']
if compiler not in compilers:
  parser.error('compiler switch "F90" has to be one out of: %s'%(', '.join(compilers)))

f2py_compiler = {
                  'gfortran': 'gnu95   --f90flags="-fPIC -fno-range-check -xf95-cpp-input -std=f2008 -fall-intrinsics -DSpectral -fdefault-real-8 -fdefault-double-8 -DFLOAT=8 -DINT=4 -I${DAMASK_ROOT}/lib"',
                  'ifort':    'intelem --f90flags="-fPIC -fpp -stand f08 -diag-disable 5268 -assume byterecl -DSpectral -real-size 64 -integer-size 32 -DFLOAT=8 -DINT=4 -I${DAMASK_ROOT}/lib"',
                }[compiler]



if   imklroot != '' and compiler != 'gfortran':
  lib_lapack = '-L%s/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -liomp5 -lpthread -lm'%(imklroot)
elif acmlroot != '':
  lib_lapack = '-L%s/%s64/lib -lacml'%(acmlroot,compiler)
elif lapackroot != '':
  lib_lapack = '-L%s -llapack'%(lapack)                                                     # see http://cens.ioc.ee/pipermail/f2py-users/2003-December/000621.html

execute = { \
          'coreModule' : [ 
                        # The following command is used to compile the fortran files and make the functions defined
                        # in damask.core.pyf available for python in the module core.so
                        # It uses the fortran wrapper f2py that is included in the numpy package to construct the
                        # module core.so out of the fortran code in the f90 files
                        # For the generation of the pyf file use the following lines:
                        ###########################################################################
                        #'f2py -h damask.core.pyf' +\
                        #' --overwrite-signature --no-lower prec.f90 DAMASK_spectral_interface.f90 math.f90 mesh.f90',
                        ###########################################################################
                        'rm `readlink -f %s`' %(os.path.join(damaskEnv.relPath('lib/damask'),'core.so')), # do this using system remove
                        'f2py damask.core.pyf' +\
                        ' --build-dir ./' +\
                        ' -c --no-lower --fcompiler=%s'%(f2py_compiler) +\
                        ' %s'%'prec.f90'+\
                        ' %s'%'DAMASK_spectral_interface.f90'+\
                        ' %s'%'IO.f90'+\
                        ' %s'%'libs.f90'+\
                        ' %s'%'numerics.f90'+\
                        ' %s'%'debug.f90'+\
                        ' %s'%'math.f90'+\
                        ' %s'%'FEsolving.f90'+\
                        ' %s'%'mesh.f90'+\
                        ' %s'%'core_quit.f90'+\
                        ' -L%s/lib -lfftw3'%(fftwroot)+\
                        ' %s'%lib_lapack,
                        'mv %s `readlink -f %s`' %(os.path.join(codeDir,'core.so'),os.path.join(damaskEnv.relPath('lib/damask'),'core.so')),    # do this using system remove
                        ]
            }

os.chdir(codeDir)                   # needed for compilation with gfortran and f2py
for tasks in execute:
  for cmd in execute[tasks]:
    try:
      print 'executing...:',cmd
      os.system(cmd)
    except:
      print 'failed..!'
      pass

modules = glob.glob('*.mod')
for module in modules:
  print 'removing', module
  os.remove(module)

#check if compilation of core module was successful
try:
  with open(damaskEnv.relPath('lib/damask/core.so')) as f: pass
except IOError as e:
  print '*********\n* core.so not found, compilation of core modules was not successful\n*********'
  sys.exit()
f.close
