#!/usr/bin/env python3
# -*- coding: UTF-8 no BOM -*-

import os
import sys
import math
import numpy as np
from optparse import OptionParser
import damask

scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

minimal_surfaces = ['primitive','gyroid','diamond',]

surface = {
            'primitive': lambda x,y,z: math.cos(x)+math.cos(y)+math.cos(z),
            'gyroid':    lambda x,y,z: math.sin(x)*math.cos(y)+math.sin(y)*math.cos(z)+math.cos(x)*math.sin(z),
            'diamond':   lambda x,y,z: math.cos(x-y)*math.cos(z)+math.sin(x+y)*math.sin(z),
          }

parser = OptionParser(option_class=damask.extendableOption, usage='%prog [option(s)] [geomfile]', description = """
Generate a geometry file of a bicontinuous structure of given type.

""", version = scriptID)


parser.add_option('-t','--type',
                  dest = 'type',
                  choices = minimal_surfaces, metavar = 'string',
                  help = 'type of minimal surface [primitive] {%s}' %(','.join(minimal_surfaces)))
parser.add_option('-f','--threshold',
                  dest = 'threshold',
                  type = 'float', metavar = 'float',
                  help = 'threshold value defining minimal surface [%default]')
parser.add_option('-g', '--grid',
                  dest = 'grid',
                  type = 'int', nargs = 3, metavar = 'int int int',
                  help = 'a,b,c grid of hexahedral box [%default]')
parser.add_option('-s', '--size',
                  dest = 'size',
                  type = 'float', nargs = 3, metavar = 'float float float',
                  help = 'x,y,z size of hexahedral box [%default]')
parser.add_option('-p', '--periods',
                  dest = 'periods',
                  type = 'int', metavar = 'int',
                  help = 'number of repetitions of unit cell [%default]')
parser.add_option('--homogenization',
                  dest = 'homogenization',
                  type = 'int', metavar = 'int',
                  help = 'homogenization index to be used [%default]')
parser.add_option('--m',
                  dest = 'microstructure',
                  type = 'int', nargs = 2, metavar = 'int int',
                  help = 'two microstructure indices to be used [%default]')
parser.set_defaults(type = minimal_surfaces[0],
                    threshold = 0.0,
                    periods = 1,
                    grid = (16,16,16),
                    size = (1.0,1.0,1.0),
                    homogenization = 1,
                    microstructure = (1,2),
                   )

(options,filenames) = parser.parse_args()

# --- loop over input files -------------------------------------------------------------------------

if filenames == []: filenames = [None]

for name in filenames:

  damask.util.report(scriptName,name)

  X = options.periods*2.0*math.pi*(np.arange(options.grid[0])+0.5)/options.grid[0]
  Y = options.periods*2.0*math.pi*(np.arange(options.grid[1])+0.5)/options.grid[1]
  Z = options.periods*2.0*math.pi*(np.arange(options.grid[2])+0.5)/options.grid[2]

  microstructure = np.empty(options.grid,dtype='int')
  for z in range(options.grid[2]):
    for y in range(options.grid[1]):
      for x in range(options.grid[0]):
        microstructure[x,y,z]=options.microstructure[options.threshold < surface[options.type](X[x],Y[y],Z[z])]
  
  geom=damask.Geom(options.size,microstructure,options.homogenization,
                   comments=[scriptID + ' ' + ' '.join(sys.argv[1:])])

  damask.util.croak('\n'.join(geom.info()))
  
  if name is None:
    sys.stdout.write(str(geom))
  else:
    geom.to_file(name)
