#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

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
""" + string.replace('$Id$','\n','\\n')
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
res=numpy.array([32,32,32])
while exitCode == 2:
  print 'restart at ', start
  proc=subprocess.Popen(executable='DAMASK_spectral',\
                        args=['-l', '%s'%options.loadcase, '-g', '%s'%options.geometry, '--regrid', '%i'%start],\
                        stderr=subprocess.PIPE,stdout=subprocess.PIPE, bufsize=1)
  while proc.poll() is None:               # while process is running
    myLine = proc.stdout.readline()
    if len(myLine)>1: print myLine[0:-1]   # print output without extra newline
  exitCode = proc.returncode
  err = proc.stderr.readlines()
  print '-------------------------------------------------------'
  print 'error messages', err
  print '-------------------------------------------------------'
  if exitCode==2:
    os.system('rm -rf %i'%start)
    os.system('mkdir %i'%start)
    os.system('cp * %i/.'%start)
    for i in xrange(len(err)):
       if re.search('restart at\s+\d+',err[i]): start=int(string.split(err[i])[2])
#------------regridding----------------------------------------------
#--------------------------------------------------------------------
    damask.core.prec.init()
    damask.core.DAMASK_interface.init(options.loadcase,options.geometry)
    damask.core.IO.init()
    damask.core.numerics.init()
    damask.core.debug.init()
    damask.core.math.init()
    damask.core.FEsolving.init()
    damask.core.mesh.init(1,1)
    damask.core.mesh.regrid(adaptive=True,resNewInput=res)
