#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,string,math
import numpy as np
from optparse import OptionParser
import damask

scriptID   = string.replace('$Id$','\n','\\n')
scriptName = os.path.splitext(scriptID.split()[1])[0]

#--------------------------------------------------------------------------------------------------
#                                MAIN
#--------------------------------------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Unpack geometry files containing ranges "a to b" and/or "n of x" multiples (exclusively in one line).
""", version = scriptID)

parser.add_option('-1', '--onedimensional',
                  dest = 'oneD',
                  action = 'store_true',
                  help = 'output geom file with one-dimensional data arrangement [%default]')

parser.set_defaults(oneD = False,
                   )

(options, filenames) = parser.parse_args()

# --- loop over input files -------------------------------------------------------------------------

if filenames == []: filenames = [None]

for name in filenames:
  try:
    table = damask.ASCIItable(name = name,
                              buffered = False, labeled = False)
  except: continue
  table.croak('\033[1m'+scriptName+'\033[0m'+(': '+name if name else ''))

# --- interpret header ----------------------------------------------------------------------------

  table.head_read()
  info,extra_header = table.head_getGeom()
  
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

# --- write header ---------------------------------------------------------------------------------

  table.labels_clear()
  table.info_clear()
  table.info_append(extra_header+[
    scriptID + ' ' + ' '.join(sys.argv[1:]),
    "grid\ta {grid[0]}\tb {grid[1]}\tc {grid[2]}".format(grid=info['grid']),
    "size\tx {size[0]}\ty {size[1]}\tz {size[2]}".format(size=info['size']),
    "origin\tx {origin[0]}\ty {origin[1]}\tz {origin[2]}".format(origin=info['origin']),
    "homogenization\t{homog}".format(homog=info['homogenization']),
    "microstructures\t{microstructures}".format(microstructures=info['microstructures']),
    ])
  table.head_write()
  
# --- write microstructure information ------------------------------------------------------------

  microstructure = table.microstructure_read(info['grid'])                                          # read microstructure
  formatwidth = int(math.floor(math.log10(microstructure.max())+1))                                 # efficient number printing format
  if options.oneD: table.data = microstructure
  else:            table.data = microstructure.reshape((info['grid'][0],info['grid'][1]*info['grid'][2]),order='F').transpose()
  table.data_writeArray('%%%ii'%(formatwidth),delimiter = ' ')
    
#--- output finalization --------------------------------------------------------------------------

  table.close()                                                                                     # close ASCII table
