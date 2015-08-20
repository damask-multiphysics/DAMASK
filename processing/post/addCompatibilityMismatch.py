#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,string
import numpy as np
from optparse import OptionParser
import damask

scriptID   = string.replace('$Id$','\n','\\n')
scriptName = os.path.splitext(scriptID.split()[1])[0]

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options file[s]', description = """
Add column(s) containing the shape and volume mismatch resulting from given deformation gradient.
Operates on periodic three-dimensional x,y,z-ordered data sets.

""", version = scriptID)


parser.add_option('-c','--coordinates',
                  dest = 'coords',
                  type = 'string', metavar = 'string',
                  help = 'column heading of coordinates [%default]')
parser.add_option('-f','--defgrad',
                  dest = 'defgrad',
                  type = 'string', metavar = 'string ',
                  help = 'column heading of deformation gradient [%default]')
parser.add_option('--no-shape','-s',
                  dest = 'shape',
                  action = 'store_false',
                  help = 'omit shape mismatch')
parser.add_option('--no-volume','-v',
                  dest = 'volume',
                  action = 'store_false',
                  help = 'omit volume mismatch')
parser.set_defaults(coords   = 'ipinitialcoord',
                    defgrad  = 'f',
                    shape = True,
                    volume = True,
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

# ------------------------------------------ assemble header --------------------------------------

  table.info_append(scriptID + '\t' + ' '.join(sys.argv[1:]))
  if options.shape:  table.labels_append('shapeMismatch({})'.format(options.defgrad))
  if options.volume: table.labels_append('volMismatch({})'.format(options.defgrad))
  table.head_write()

# --------------- figure out size and grid ---------------------------------------------------------

  table.data_readArray()

  coords = [np.unique(table.data[:,colCoord+i]) for i in xrange(3)]
  mincorner = np.array(map(min,coords))
  maxcorner = np.array(map(max,coords))
  grid   = np.array(map(len,coords),'i')
  size   = grid/np.maximum(np.ones(3,'d'), grid-1.0) * (maxcorner-mincorner)                        # size from edge to edge = dim * n/(n-1) 
  size   = np.where(grid > 1, size, min(size[grid > 1]/grid[grid > 1]))                             # spacing for grid==1 equal to smallest among other spacings

  N = grid.prod()
  
# ------------------------------------------ process deformation gradient --------------------------

  F = table.data[:,colF:colF+9].transpose().reshape([3,3]+list(options.dimension),order='F')
  Favg    = damask.core.math.tensorAvg(F)
  centres = damask.core.mesh.deformedCoordsFFT(size,F,Favg,[1.0,1.0,1.0])
  nodes   = damask.core.mesh.nodesAroundCentres(size,Favg,centres)

  stack = [table.data]
  if options.shape:  stack.append(damask.core.mesh.shapeMismatch( size,F,nodes,centres))
  if options.volume: stack.append(damask.core.mesh.volumeMismatch(size,F,nodes))

# ------------------------------------------ output result -----------------------------------------

  if len(stack) > 1: table.data = np.hstack(tuple(stack))
  table.data_writeArray()

# ------------------------------------------ output finalization -----------------------------------  

  table.close()                                                                                     # close ASCII tables
