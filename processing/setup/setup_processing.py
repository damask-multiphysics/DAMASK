#!/usr/bin/env python

# Makes postprocessing routines acessible from everywhere.

import os,sys,glob,string
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

""" + string.replace('$Id: addDivergence.py 1129 2011-12-01 12:01:13Z MPIE\m.diehl $','\n','\\n')
)

compilers = ['intel','ifort','intel32','gfortran','gnu95']

parser.add_option('--F90',              dest='compiler', type='string', \
                                        help='name of F90 compiler')
                                        
(options,filenames) = parser.parse_args()

if options.compiler not in compilers:
  parser.error('compiler switch "--F90" has to be one out of: %s'%(', '.join(compilers)))

f2py_compiler = {
                  'gfortran': 'gnu95  --f90flags="-fno-range-check"',
                  'gnu95':    'gnu95  --f90flags="-fno-range-check"',
                  'intel32':  'intel',
                  'intel':    'intelem',
                  'ifort':    'intelem',
                }[options.compiler]
compiler = {
                  'gfortran': 'gfortran',
                  'gnu95':    'gfortran',
                  'intel32':  'ifort',
                  'intel':    'ifort',
                  'ifort':    'ifort',
                }[options.compiler]

damaskEnv = Environment()
baseDir = damaskEnv.relPath('processing/')
codeDir = damaskEnv.relPath('code/')

if   'ikml' in damaskEnv.pathInfo and damaskEnv.pathInfo['ikml'] != '':
  lib_lapack = ''  # TODO!!
elif 'acml' in damaskEnv.pathInfo and damaskEnv.pathInfo['acml'] != '':
  lib_lapack = '-L%s -lacml'%(os.path.join(damaskEnv.pathInfo['acml'],'%s64'))                                    # can we use linker flag?
#  lib_lapack = os.path.join(damaskEnv.pathInfo['acml'],'%s64/lib/libacml.a'%options.compiler)     # why linking against static lib?
elif 'lapack' in damaskEnv.pathInfo and damaskEnv.pathInfo['lapack'] != '':
  lib_lapack = '-L%s -llapack'%(damaskEnv.pathInfo['lapack'])                                     # see http://cens.ioc.ee/pipermail/f2py-users/2003-December/000621.html

#define ToDo list
bin_link = { \
        'pre' : [
                'marc_addUserOutput.py',
                'mentat_pbcOnBoxMesh.py',
                'mentat_spectralBox.py',
                'patchFromReconstructedBoundaries.py',
                'spectral_geomCheck.py',
                'spectral_geomScale.py',
                'spectral_geomCrop.py',
                'spectral_minimalSurface.py',
                'spectral_vicinityOffset.py',
                'voronoi_randomSeeding.exe',
                'voronoi_tessellation.exe',
                ],
        'post' : [
                '3Dvisualize.py',
                'addCauchy.py',
                'addCalculation.py',
                'addDeterminant.py',
                'addDivergence.py',
                'addCurl.py',
                'addMises.py',
                'addNorm.py',
                'addStrainTensors.py',
                'addCompatibilityMismatch.py',
                'averageDown.py',
                'mentat_colorMap.py',
                'postResults.py',
                'spectral_iterationCount.py',
                'spectral_convergence.py',
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
          'postMath' : [ 
                        'make tidy',
                        'rm %s'%(os.path.join(damaskEnv.relPath('lib/damask'),'core.so')),
                        # The following command is used to compile math.f90 and make the functions defined in DAMASK_math.pyf
                        # available for python in the module DAMASK_math.so
                        # It uses the fortran wrapper f2py that is included in the numpy package to construct the
                        # module postprocessingMath.so out of the fortran code postprocessingMath.f90
                        # for the generation of the pyf file:
                        #f2py -m DAMASK -h DAMASK.pyf --overwrite-signature ../../code/math.f90 \
                        'f2py %s'%(os.path.join(codeDir,'damask.core.pyf')) +\
                        ' -c --fcompiler=%s'%(f2py_compiler) +\
                        ' %s'%(os.path.join(codeDir,'core_modules.f90'))+\
                        ' %s'%(os.path.join(codeDir,'math.f90'))+\
                        ' -L%s -lfftw3'%(damaskEnv.pathInfo['fftw'])+\
                        ' %s'%lib_lapack,
                        'mv %s %s' %(os.path.join(codeDir,'core.so'),damaskEnv.relPath('lib/damask')),
                        ]
            }


for dir in compile:
  for file in compile[dir]:
    src = os.path.abspath(os.path.join(baseDir,dir,file))
    if os.path.isfile(src):
      try:
        os.system('rm %s.exe'%(os.path.splitext(src)[0]))
        print 'removing %s.exe '%(os.path.splitext(src)[0])
      except:
        pass
      print 'compiling ',src,'using',compiler
      os.system('%s -O2 -o%s.exe %s'%(compiler,os.path.splitext(src)[0],src))
 
os.chdir(codeDir)                   # needed for compilation with gfortran and f2py
for tasks in execute:
  for cmd in execute[tasks]:
    try:
      print 'executing...:',cmd
      os.system(cmd)
    except:
      print 'failed..!'
      pass

os.chdir(damaskEnv.relPath('processing/setup/'))
modules = glob.glob('*.mod')
for module in modules:
  print 'removing', module
  os.remove(module)
  
for dir in bin_link:
  for file in bin_link[dir]:
    src = os.path.abspath(os.path.join(baseDir,dir,file))
    if (file == ''):
      sym_link = os.path.abspath(os.path.join(damaskEnv.binDir(),dir))
    else:
      sym_link = os.path.abspath(os.path.join(damaskEnv.binDir(),os.path.splitext(file)[0]))
    print sym_link,'-->',src
    if os.path.lexists(sym_link):
      os.remove(sym_link)    
    os.symlink(src,sym_link)
    
