#!/usr/bin/env python

# compiles fortran code for Python

import os,sys,glob,string,subprocess,shlex
from damask import Environment

damaskEnv = Environment()
baseDir = damaskEnv.relPath('installation/')
codeDir = damaskEnv.relPath('code/')

options={}

for option in ['IMKL_ROOT','ACML_ROOT','LAPACK_ROOT','FFTW_ROOT','F90']:
  try:
    value = damaskEnv.options[option]
  except:
    value = os.getenv(option)
    if value is None: value = ''           # env not set
  options[option]=value

compilers = ['ifort','gfortran']
if options['F90'] not in compilers:
  sys.exit('compiler "F90" (in installation/options or as Shell variable) has to be one out of: %s'%(', '.join(compilers)))
print 'Argument List:', str(sys.argv)
compileCommand = {
                  'gfortran': 'gnu95 --f90flags="-fPIC -fno-range-check -xf95-cpp-input -std=f2008 -fall-intrinsics'+\
                          ' -DSpectral -fdefault-real-8 -fdefault-double-8 -DFLOAT=8 -DINT=4 -I%s/lib"'%damaskEnv.rootDir(),
                  'ifort':    'intelem --f90flags="-fPIC -fpp -stand f08 -diag-disable 5268 -assume byterecl'+\
                          ' -DSpectral -real-size 64 -integer-size 32 -DFLOAT=8 -DINT=4 -I%s/lib"'%damaskEnv.rootDir(),
                  }[options['F90']]

# see http://cens.ioc.ee/pipermail/f2py-users/2003-December/000621.html
if   options['IMKL_ROOT'] != '' and options['F90'] != 'gfortran':
  lib_lapack = '-L%s/lib/intel64 -I%s/include -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm -liomp5'%(options['IMKL_ROOT'],options['IMKL_ROOT'])
elif options['ACML_ROOT'] != '':
  lib_lapack = '-L%s/%s64/lib -lacml'%(options['ACML_ROOT'],options['F90'])
elif options['LAPACK_ROOT'] != '':
  lib_lapack = '-L%s/lib -L%s/lib64 -llapack'%(options['LAPACK_ROOT'],options['LAPACK_ROOT'])                                                     

os.chdir(codeDir)                                                                           # needed for compilation with gfortran and f2py
try:
  os.remove(os.path.join(damaskEnv.relPath('lib/damask'),'core.so'))
except OSError, e:                                                                          ## if failed, report it back to the user ##
  print ("Error when deleting: %s - %s." % (e.filename,e.strerror))


# The following command is used to compile the fortran files and make the functions defined
# in damask.core.pyf available for python in the module core.so
# It uses the fortran wrapper f2py that is included in the numpy package to construct the
# module core.so out of the fortran code in the f90 files
# For the generation of the pyf file use the following lines:
###########################################################################
#'f2py -h damask.core.pyf' +\
#' --overwrite-signature --no-lower prec.f90 DAMASK_spectral_interface.f90 math.f90 mesh.f90,...'
 ###########################################################################
cmd = 'f2py damask.core.pyf' +\
      ' -c --no-lower --fcompiler=%s'%(compileCommand) +\
      ' prec.f90'+\
      ' DAMASK_spectral_interface.f90'+\
      ' IO.f90'+\
      ' libs.f90'+\
      ' numerics.f90'+\
      ' debug.f90'+\
      ' math.f90'+\
      ' FEsolving.f90'+\
      ' mesh.f90'+\
      ' core_quit.f90'+\
      ' -L%s/lib -lfftw3'%(options['FFTW_ROOT'])+\
      ' %s'%lib_lapack

print('Executing: '+cmd)
try:
  subprocess.call(shlex.split(cmd))
except subprocess.CalledProcessError:
  print('build failed')
except OSError:
  print ('f2py not found')

try:
  os.rename(os.path.join(codeDir,'core.so'),\
            os.path.join(damaskEnv.relPath('lib/damask'),'core.so'))
except:
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
