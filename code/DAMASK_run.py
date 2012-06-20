#!/usr/bin/env python

import numpy,os,damask,string,sys,subprocess,re
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

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=extendableOption, usage='%prog options [file[s]]', description = """
Add column(s) with derived values according to user defined arithmetic operation between column(s).
Columns can be specified either by label or index.

Example: distance to IP coordinates -- "math.sqrt( #ip.x#**2 + #ip.y#**2 + #ip.z#**2 )"
""" + string.replace('$Id: addCalculation.py 1355 2012-02-23 13:53:12Z MPIE\p.eisenlohr $','\n','\\n')
)

parser.add_option('-l','--load', '--loadcase',   dest='loadcase',  type='string', \
                                    help='PathToLoadFile/NameOfLoadFile.load. "PathToLoadFile" will be the working directory.')
parser.add_option('-g','--geom', '--geometry', dest='geometry', type='string', \
                                    help='PathToGeomFile/NameOfGeomFile.load.')

parser.set_defaults(loadcase= '')
parser.set_defaults(geometry= '')

(options,filenames) = parser.parse_args()
start = 1
exitCode=2
print 'load case', options.loadcase
print 'geometry', options.geometry
f=open('monitor','w')
#res=numpy.array([6,6,6])
res=numpy.array([0,0,0])
while exitCode == 2:
  print 'restart at ', start
  out=subprocess.Popen(['DAMASK_spectral', '-l', '%s'%options.loadcase, '-g', '%s'%options.geometry, '--regrid', '%i'%start],stderr=subprocess.PIPE,stdout=f)
  #out=subprocess.Popen(['DAMASK_spectral', '-l', '%s'%options.loadcase, '-g', '%s'%options.geometry, '-r', '%i'%start],stderr=subprocess.PIPE,stdout=f)
  stderr = out.communicate()
  stderrLines = string.split(stderr[1],'\n')
  exitCode = int(stderrLines[-2])
  print 'exit code', exitCode
  if exitCode==2:
    os.system('rm -rf %i'%start)
    os.system('mkdir %i'%start)
    os.system('cp * %i/.'%start)
    start = int(string.split(re.search('restart at\s+\d+',stderr[1]).group(0))[2])
#------------regridding----------------------------------------------
#--------------------------------------------------------------------
    damask.core.prec.prec_init()
    damask.core.damask_interface.damask_interface_init(options.loadcase,options.geometry)
    damask.core.io.io_init()
    damask.core.numerics.numerics_init()
    damask.core.debug.debug_init()
    damask.core.math.math_init()
    damask.core.fesolving.fe_init()
    damask.core.mesh.mesh_init(1,1)
    damask.core.mesh.mesh_regrid(resNewInput=res)
