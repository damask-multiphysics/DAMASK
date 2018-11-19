#!/usr/bin/env python3
# -*- coding: UTF-8 no BOM -*-

import os,sys,math
import numpy as np
import damask
from optparse import OptionParser

scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])

#--------------------------------------------------------------------------------------------------
#                                MAIN
#--------------------------------------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog [file[s]]', description = """
renumber sorted microstructure indices to 1,...,N.

""", version=scriptID)

(options, filenames) = parser.parse_args()

# --- loop over input files ----------------------------------------------------------------------

if filenames == []: filenames = [None]

for name in filenames:
  try:    table = damask.ASCIItable(name = name,
                                    buffered = False,
                                    labeled = False)
  except: continue
  damask.util.report(scriptName,name)

# --- interpret header ---------------------------------------------------------------------------

  table.head_read()
  info,extra_header = table.head_getGeom()

  damask.util.croak(['grid     a b c:  %s'%(' x '.join(map(str,info['grid']))),
               'size     x y z:  %s'%(' x '.join(map(str,info['size']))),
               'origin   x y z:  %s'%(' : '.join(map(str,info['origin']))),
               'homogenization:  %i'%info['homogenization'],
               'microstructures: %i'%info['microstructures'],
              ])

  errors = []
  if np.any(info['grid'] < 1):    errors.append('invalid grid a b c.')
  if np.any(info['size'] <= 0.0): errors.append('invalid size x y z.')
  if errors != []:
    damask.util.croak(errors)
    table.close(dismiss = True)
    continue

# --- read data ----------------------------------------------------------------------------------

  microstructure = table.microstructure_read(info['grid'])                                 # read microstructure

# --- do work ------------------------------------------------------------------------------------

  newInfo = {
             'origin':  np.zeros(3,'d'),
             'microstructures': 0,
            }

  grainIDs = np.unique(microstructure)
  renumbered = np.copy(microstructure)
  
  for i, oldID in enumerate(grainIDs):
    renumbered = np.where(microstructure == oldID, i+1, renumbered)

  newInfo['microstructures'] = len(grainIDs)

# --- report -------------------------------------------------------------------------------------

  remarks = []
  if (    newInfo['microstructures'] != info['microstructures']):
    remarks.append('--> microstructures: %i'%newInfo['microstructures'])
  if remarks != []: damask.util.croak(remarks)

# --- write header -------------------------------------------------------------------------------

  table.labels_clear()
  table.info_clear()
  table.info_append(extra_header+[
    scriptID + ' ' + ' '.join(sys.argv[1:]),
    "grid\ta {grid[0]}\tb {grid[1]}\tc {grid[2]}".format(grid=info['grid']),
    "size\tx {size[0]}\ty {size[1]}\tz {size[2]}".format(size=info['size']),
    "origin\tx {origin[0]}\ty {origin[1]}\tz {origin[2]}".format(origin=info['origin']),
    "homogenization\t{homog}".format(homog=info['homogenization']),
    "microstructures\t{microstructures}".format(microstructures=newInfo['microstructures']),
    ])
  table.head_write()

# --- write microstructure information -----------------------------------------------------------

  format = '%{}i'.format(int(math.floor(math.log10(newInfo['microstructures'])+1)))
  table.data = renumbered.reshape((info['grid'][0],info['grid'][1]*info['grid'][2]),order='F').transpose()
  table.data_writeArray(format,delimiter = ' ')

# --- output finalization ------------------------------------------------------------------------

  table.close()                                                                                   # close ASCII table
