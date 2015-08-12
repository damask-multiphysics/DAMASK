#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,string,math
import numpy as np
from optparse import OptionParser
import damask

scriptID   = string.replace('$Id$','\n','\\n')
scriptName = os.path.splitext(scriptID.split()[1])[0]

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

minimal_surfaces = ['primitive','gyroid','diamond',]

surface = {
            'primitive': lambda x,y,z: math.cos(x)+math.cos(y)+math.cos(z),
            'gyroid':    lambda x,y,z: math.sin(x)*math.cos(y)+math.sin(y)*math.cos(z)+math.cos(x)*math.sin(z),
            'diamond':   lambda x,y,z: math.cos(x-y)*math.cos(z)+math.sin(x+y)*math.sin(z),
          }

parser = OptionParser(option_class=damask.extendableOption, usage='%prog', description = """
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
parser.add_option('-1', '--onedimensional',
                  dest = 'oneD',
                  action = 'store_true',
                  help = 'output geom file with two-dimensional data arrangement [%default]')
parser.set_defaults(type = minimal_surfaces[0],
                    threshold = 0.0,
                    periods = 1,
                    grid = (16,16,16),
                    size = (1.0,1.0,1.0),
                    homogenization = 1,
                    microstructure = (1,2),
                    oneD  = False,
                   )

(options,filenames) = parser.parse_args()

# --- loop over input files -------------------------------------------------------------------------

if filenames == []: filenames = [None]

for name in filenames:
  try:
    table = damask.ASCIItable(outname = name,
                              buffered = False, labeled = False)
  except:
    continue
  table.croak('\033[1m'+scriptName+'\033[0m'+(': '+name if name else ''))


# ------------------------------------------ make grid -------------------------------------

  info = {
          'grid':   np.array(options.grid),
          'size':   np.array(options.size),
          'origin': np.zeros(3,'d'),
          'microstructures': max(options.microstructure),
          'homogenization':  options.homogenization
         }

#--- report ---------------------------------------------------------------------------------------

  table.croak(['grid     a b c:  %s'%(' x '.join(map(str,info['grid']))),
               'size     x y z:  %s'%(' x '.join(map(str,info['size']))),
               'origin   x y z:  %s'%(' : '.join(map(str,info['origin']))),
               'homogenization:  %i'%info['homogenization'],
               'microstructures: %i'%info['microstructures'],
              ])

  errors = []
  if np.any(info['grid'] < 1):    errors.append('invalid grid a b c.')
  if np.any(info['size'] <= 0.0): errors.append('invalid size x y z.')
  if errors != []:
    table.croak(errors)
    table.close(dismiss = True)
    continue

#--- write header ---------------------------------------------------------------------------------

  table.labels_clear()
  table.info_clear()
  table.info_append([
    scriptID + ' ' + ' '.join(sys.argv[1:]),
    "grid\ta {grid[0]}\tb {grid[1]}\tc {grid[2]}".format(grid=info['grid']),
    "size\tx {size[0]}\ty {size[1]}\tz {size[2]}".format(size=info['size']),
    "origin\tx {origin[0]}\ty {origin[1]}\tz {origin[2]}".format(origin=info['origin']),
    "homogenization\t{homog}".format(homog=info['homogenization']),
    "microstructures\t{microstructures}".format(microstructures=info['microstructures']),
    ])
  table.head_write()
  
#--- write data -----------------------------------------------------------------------------------
  X = options.periods*2.0*math.pi*(np.arange(options.grid[0])+0.5)/options.grid[0]
  Y = options.periods*2.0*math.pi*(np.arange(options.grid[1])+0.5)/options.grid[1]
  Z = options.periods*2.0*math.pi*(np.arange(options.grid[2])+0.5)/options.grid[2]

  for z in xrange(options.grid[2]):
    for y in xrange(options.grid[1]):
      table.data_clear()
      for x in xrange(options.grid[0]):
        table.data_append(options.microstructure[options.threshold < surface[options.type](X[x],Y[y],Z[z])])
        if options.oneD:
          table.data_write()
          table.data_clear()
      table.data_write()

  table.close()
