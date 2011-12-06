#!/usr/bin/env python

# Makes postprocessing routines acessible from everywhere.

import os,sys,glob,string
from optparse import OptionParser, Option

import damask_tools

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

parser.add_option('--F90',              dest='compiler', type='string', \
                                        help='name of F90 compiler [%default]')
                                        
parser.set_defaults(compiler = 'ifort')
(options,filenames) = parser.parse_args()

#translating name of compiler for use with f2py and setting subdirname of acml
if options.compiler == 'gfortran':
  f2py_compiler='gnu95  --f90flags="-fno-range-check"'
else:
  f2py_compiler='intelem'

acml_subdir='%s64/lib'%options.compiler


#getting pathinfo
damask_variables = damask_tools.DAMASK_TOOLS()
baseDir = damask_variables.relPath('processing/')
codeDir = damask_variables.relPath('code/')

#define ToDo list
bin_link = { \
        'pre' : [
                'marc_addUserOutput.py',
                'mentat_pbcOnBoxMesh.py',
                'mentat_spectralBox.py',
                'patchFromReconstructedBoundaries.py',
                'spectral_geomCheck.py',
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
          'postMath' : [ 
                        'rm %s'%(os.path.join(damask_variables.relPath('lib/'),'DAMASK.so')),
                        # The following command is used to compile math.f90 and make the functions defined in DAMASK_math.pyf
                        # available for python in the module DAMASK_math.so
                        # It uses the fortran wrapper f2py that is included in the numpy package to construct the
                        # module postprocessingMath.so out of the fortran code postprocessingMath.f90
                        # for the generation of the pyf file:
                        #f2py -m DAMASK -h DAMASK.pyf --overwrite-signature ../../code/math.f90 \
                        'f2py %s '%(os.path.join(codeDir,'DAMASK.pyf')) +\
                        '-c --fcompiler=%s '%(f2py_compiler) +\
                        '%s ' %(os.path.join(codeDir,'DAMASK2Python_helper.f90'))+\
                        '%s ' %(os.path.join(codeDir,'math.f90'))+\
                        '%s ' %(os.path.join(damask_variables.pathInfo['fftw'],'libfftw3.a'))+\
                        '%s' %(os.path.join(damask_variables.pathInfo['acml'],acml_subdir,'libacml.a'),
                        'mv %s %s' %(os.path.join(codeDir,'DAMASK.so'),damask_variables.relPath('lib/')),
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
      print 'compiling ',src,'using',options.compiler
      os.system('%s -O2 -o%s.exe %s'%(options.compiler,os.path.splitext(src)[0],src))
 
os.chdir(codeDir)                   # needed for compilation with gfortran and f2py
for tasks in execute:
  for cmd in execute[tasks]:
    os.system(cmd)  
os.chdir(damask_variables.relPath('processing/setup/'))

modules = glob.glob('*.mod')
for module in modules:
  print 'removing', module
  os.remove(module)
  
for dir in bin_link:
  for file in bin_link[dir]:
    src = os.path.abspath(os.path.join(baseDir,dir,file))
    if (file == ''):
      sym_link = os.path.abspath(os.path.join(damask_variables.binDir(),dir))
    else:
      sym_link = os.path.abspath(os.path.join(damask_variables.binDir(),os.path.splitext(file)[0]))
    print sym_link,'-->',src
    if os.path.lexists(sym_link):
      os.remove(sym_link)    
    os.symlink(src,sym_link)
    
    #--- uncomment next lines to remove your old symbolic links in ~/bin
    #old_link=sym_link.replace(damask_root,os.getenv('HOME'))
    #if os.path.lexists(old_link):
    #  os.remove(old_link)

