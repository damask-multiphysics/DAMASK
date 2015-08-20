#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,string
import numpy as np
from collections import defaultdict
from optparse import OptionParser
import damask

scriptID   = string.replace('$Id$','\n','\\n')
scriptName = os.path.splitext(scriptID.split()[1])[0]

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options file[s]', description = """
Add deformed configuration of given initial coordinates.
Operates on periodic three-dimensional x,y,z-ordered data sets.

""", version = scriptID)

parser.add_option('-c','--coordinates',
                  dest = 'coords',
                  type = 'string', metavar = 'string',
                  help = 'column label of coordinates [%default]')
parser.add_option('-f','--defgrad',
                  dest = 'defgrad',
                  type = 'string', metavar = 'string',
                  help = 'column label of deformation gradient [%default]')
parser.add_option('--scaling',
                  dest = 'scaling',
                  type = 'float', nargs = 3, metavar = ' '.join(['float']*3),
                  help = 'x/y/z scaling of displacement fluctuation')

parser.set_defaults(coords  = 'ipinitialcoord',
                    defgrad = 'f',
                    scaling = [1.,1.,1.],
                   )

(options,filenames) = parser.parse_args()

# --- loop over input files -------------------------------------------------------------------------

if filenames == []: filenames = [None]

for name in filenames:
  try:
    table = damask.ASCIItable(name = name,
                              buffered = False)
  except: continue
  table.croak('\033[1m'+scriptName+'\033[0m'+(': '+name if name else ''))

# ------------------------------------------ read header ------------------------------------------

  table.head_read()

# ------------------------------------------ sanity checks ----------------------------------------

  errors  = []
  remarks = []
  
  if table.label_dimension(options.coords) != 3:  errors.append('coordinates {} are not a vector.'.format(options.coords))
  else: colCoord = table.label_index(options.coords)

  if table.label_dimension(options.defgrad) != 9: errors.append('deformation gradient {} is not a tensor.'.format(options.defgrad))
  else: colF = table.label_index(options.defgrad)

  if remarks != []: table.croak(remarks)
  if errors  != []:
    table.croak(errors)
    table.close(dismiss = True)
    continue

# --------------- figure out size and grid ---------------------------------------------------------

  table.data_readArray()

  coords = [np.unique(table.data[:,colCoord+i]) for i in xrange(3)]
  mincorner = np.array(map(min,coords))
  maxcorner = np.array(map(max,coords))
  grid   = np.array(map(len,coords),'i')
  size   = grid/np.maximum(np.ones(3,'d'), grid-1.0) * (maxcorner-mincorner)                        # size from edge to edge = dim * n/(n-1) 
  size   = np.where(grid > 1, size, min(size[grid > 1]/grid[grid > 1]))                             # spacing for grid==1 equal to smallest among other spacings

  N = grid.prod()
  
  if N != len(table.data): errors.append('data count {} does not match grid {}x{}x{}.'.format(N,*grid))
  if errors  != []:
    table.croak(errors)
    table.close(dismiss = True)
    continue
  
# ------------------------------------------ assemble header ---------------------------------------

  table.info_append(scriptID + '\t' + ' '.join(sys.argv[1:]))
  table.labels_append(['{}_{}.{}'%(coord+1,options.defgrad,options.coords) for coord in xrange(3)])  # extend ASCII header with new labels
  table.head_write()

# ------------------------------------------ process deformation gradient --------------------------

  F = table.data[:,colF:colF+9].transpose().reshape([3,3]+list(options.dimension),order='F')
  Favg    = damask.core.math.tensorAvg(F)
  centres = damask.core.mesh.deformedCoordsFFT(size,F,Favg,[1.0,1.0,1.0])

  stack = [table.data,centres]

# ------------------------------------------ output result -----------------------------------------

  if len(stack) > 1: table.data = np.hstack(tuple(stack))
  table.data_writeArray()

# ------------------------------------------ output finalization -----------------------------------  

  table.close()                                                                                     # close ASCII tables
