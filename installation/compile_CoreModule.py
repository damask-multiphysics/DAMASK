#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,glob,string,subprocess,shlex
from damask import Environment
from damask import version as DAMASKVERSION

# compiles fortran code for Python
scriptID = '$Id$'

damaskEnv = Environment()
baseDir = damaskEnv.relPath('installation/')
codeDir = damaskEnv.relPath('code/')

keywords=['IMKL_ROOT','ACML_ROOT','LAPACK_ROOT','FFTW_ROOT','F90']
options={}

#--- getting options from damask.conf or, if not present, from envinronment -----------------------
for option in keywords:
  try:
    value = damaskEnv.options[option]
  except:
    value = os.getenv(option)
    if value is None: value = ''           # env not set
  options[option]=value

#--- overwrite default options with keyword=value pair from argument list to mimic make behavior --
for i, arg in enumerate(sys.argv):
  for option in keywords:
    if arg.startswith(option):
      options[option] = sys.argv[i][len(option)+1:]

#--- check for valid compiler and set options -----------------------------------------------------
compilers = ['ifort','gfortran']
if options['F90'] not in compilers:
  sys.exit('compiler "F90" (in installation/options or as Shell variable) has to be one out of: %s'%(', '.join(compilers)))

compiler = {
            'gfortran': '--fcompiler=gnu95 --f90flags="-fPIC -fno-range-check -xf95-cpp-input -std=f2008 -fall-intrinsics'+\
                        ' -fdefault-real-8 -fdefault-double-8"',
            'ifort':    '--fcompiler=intelem --f90flags="-fPIC -fpp -stand f08 -diag-disable 5268 -assume byterecl'+\
                        ' -real-size 64 -integer-size 32 -shared-intel"',
            }[options['F90']]

#--- option not depending on compiler -------------------------------------------------------------
compileOptions = ' -DSpectral -DFLOAT=8 -DINT=4 -I%s/lib -DDAMASKVERSION=\\\\\"\\\"%s\\\\\"\\\"'%(damaskEnv.rootDir(),DAMASKVERSION)

#--- this saves the path of libraries to core.so, hence it is known during runtime ----------------
if options['F90'] == 'gfortran':
  LDFLAGS = '-shared -Wl,-undefined,dynamic_lookup'                                                 # solved error: Undefined symbols for architecture x86_64: "_PyArg_ParseTupleAndKeywords" as found on https://lists.macosforge.org/pipermail/macports-dev/2013-May/022735.html
else:
  LDFLAGS = ' -openmp -Wl'                                                                          # some f2py versions/configurations compile with openMP, so linking against openMP is needed to prevent errors during loading of core module

#--- run path of for fftw during runtime ----------------------------------------------------------
LDFLAGS += ',-rpath,%s/lib,-rpath,%s/lib64'%(options['FFTW_ROOT'],options['FFTW_ROOT'])

# see http://cens.ioc.ee/pipermail/f2py-users/2003-December/000621.html
if options['IMKL_ROOT']:
  if options['F90'] == 'gfortran':
    arch = 'gf'
  elif options['F90'] == 'ifort':
    arch = 'intel'
  lib_lapack = '-L%s/lib/intel64 -lmkl_%s_lp64 -lmkl_core -lmkl_sequential -lm'\
                                                                      %(options['IMKL_ROOT'],arch)
  LDFLAGS +=',-rpath,%s/lib/intel64'%(options['IMKL_ROOT'])
elif options['ACML_ROOT'] != '':
  lib_lapack = '-L%s/%s64/lib  -lacml'%(options['ACML_ROOT'],options['F90'])
  LDFLAGS +=',-rpath,%s/%s64/lib'%(options['ACML_ROOT'],options['F90'])
elif options['LAPACK_ROOT'] != '':
  lib_lapack = '-L%s/lib -L%s/lib64  -llapack'%(options['LAPACK_ROOT'],options['LAPACK_ROOT'])
  LDFLAGS +=',-rpath,%s/lib,-rpath,%s/lib64'%(options['LAPACK_ROOT'],options['LAPACK_ROOT'])

#--------------------------------------------------------------------------------------------------
# f2py does not (yet) support setting of special flags for the linker, hence they must be set via 
# environment variable ----------------------------------------------------------------------------
my_env = os.environ
my_env["LDFLAGS"] = LDFLAGS                                    

#--- delete old file ------------------------------------------------------------------------------
try:
  os.remove(os.path.join(damaskEnv.relPath('lib/damask'),'core.so'))
except OSError, e:
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
os.chdir(codeDir)
cmd = 'f2py damask.core.pyf' +\
      ' -c --no-lower %s'%(compiler) +\
      compileOptions+\
      ' prec.f90'+\
      ' spectral_interface.f90'+\
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
  subprocess.call(shlex.split(cmd),env=my_env)
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

#--- check if compilation of core module was successful -------------------------------------------
try:
  with open(damaskEnv.relPath('lib/damask/core.so')) as f: pass
except IOError as e:
  print '*********\n* core.so not found, compilation of core modules was not successful\n*********'
