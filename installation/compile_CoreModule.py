#!/usr/bin/env python

# Makes postprocessing routines acessible from everywhere.

import os,sys,glob,string,subprocess
from damask import Environment
from optparse import OptionParser, Option


# -----------------------------
class extendableOption(Option):
# -----------------------------
# used for definition of new option parser action 'extend', which enables to take multiple option arguments
# taken from online tutorial http://docs.python.org/library/optparse.html
  
  ACTIONS = Option.ACTIONS + ("extend",)
  STORE_ACTIONS = Option.STORE_ACTIONS + ("extend",)
  TYPED_ACTIONS = Option.TYPED_ACTIONS + ("extend",)
  ALWAYS_TYPED_ACTIONS = Option.ALWAYS_TYPED_ACTIONS + ("extend",)

  def take_action(self, action, dest, opt, value, values, parser):
    if action == "extend":
      lvalue = value.split(",")
      values.ensure_value(dest, []).extend(lvalue)
    else:
      Option.take_action(self, action, dest, opt, value, values, parser)

      
      
########################################################
# MAIN
########################################################

parser = OptionParser(option_class=extendableOption, usage='%prog options', description = """
Sets up the pre and post processing tools of DAMASK

""" + string.replace('$Id$','\n','\\n')
)

compilers = ['ifort','gfortran']

parser.add_option('--F90',          dest='compiler', type='string', \
                                    help='name of F90 [%default]')
parser.add_option('--FFTWROOT',     dest='fftwroot', type='string', \
                                    help='root location of FFTW [%default]')
parser.add_option('--LAPACKROOT',   dest='lapackroot', type='string', \
                                    help='root location of LAPACK [%default]')
parser.add_option('--ACMLROOT',     dest='acmlroot', type='string', \
                                    help='root location of ACML [%default]')
parser.add_option('--IMKLROOT',     dest='imklroot', type='string', \
                                    help='root location of IMKL [%default]')

parser.set_defaults(compiler   = {True:os.getenv('F90'),False:''}\
                                  [os.getenv('F90')!=None])
parser.set_defaults(fftwroot   = {True:os.getenv('FFTWROOT'),False:''}\
                                  [os.getenv('FFTWROOT')!=None])
parser.set_defaults(lapackroot = {True:os.getenv('LAPACKROOT'),False:''}\
                                  [os.getenv('LAPACKROOT')!=None])
parser.set_defaults(acmlroot   = {True:os.getenv('ACMLROOT'),False:''}\
                                  [os.getenv('ACMLROOT')!=None])
parser.set_defaults(imklroot   = {True:os.getenv('IMKLROOT'),False:''}\
                                  [os.getenv('IMKLROOT')!=None])

(options,filenames) = parser.parse_args()

if options.compiler not in compilers:
  parser.error('compiler switch "F90" has to be one out of: %s'%(', '.join(compilers)))

f2py_compiler = {
                  'gfortran': 'gnu95   --f90flags="-fPIC -fno-range-check -xf95-cpp-input -std=f2008 -fall-intrinsics -DSpectral -fdefault-real-8 -fdefault-double-8 -DFLOAT=8 -DINT=4 -I${DAMASK_ROOT}/lib"',
                  'ifort':    'intelem --f90flags="-fPIC -fpp -stand f08 -diag-disable 5268 -assume byterecl -DSpectral -real-size 64 -integer-size 32 -DFLOAT=8 -DINT=4 -I${DAMASK_ROOT}/lib"',
                }[options.compiler]

damaskEnv = Environment()
baseDir = damaskEnv.relPath('installation/')
codeDir = damaskEnv.relPath('code/')

if   options.imklroot != '' and options.compiler != 'gfortran':
  lib_lapack = '-L%s/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -liomp5 -lpthread -lm'%(options.imklroot)
elif options.acmlroot != '':
  lib_lapack = '-L%s/%s64/lib -lacml'%(options.acmlroot,options.compiler)
elif options.lapackroot != '':
  lib_lapack = '-L%s -llapack'%(options.lapack)                                                     # see http://cens.ioc.ee/pipermail/f2py-users/2003-December/000621.html

execute = { \
          'coreModule' : [ 
                        'make tidy',
                        # The following command is used to compile the fortran files and make the functions defined
                        # in damask.core.pyf available for python in the module core.so
                        # It uses the fortran wrapper f2py that is included in the numpy package to construct the
                        # module core.so out of the fortran code in the f90 files
                        # For the generation of the pyf file use the following lines:
                        ###########################################################################
                        #'f2py -h damask.core.pyf' +\
                        #' --overwrite-signature --no-lower prec.f90 DAMASK_spectral_interface.f90 math.f90 mesh.f90',
                        ###########################################################################
                        'rm `readlink -f %s`' %(os.path.join(damaskEnv.relPath('lib/damask'),'core.so')),
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
                        ' -L%s/lib -lfftw3'%(options.fftwroot)+\
                        ' %s'%lib_lapack,
                        'mv %s `readlink -f %s`' %(os.path.join(codeDir,'core.so'),os.path.join(damaskEnv.relPath('lib/damask'),'core.so')),
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
