#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,string
import numpy as np
import scipy.ndimage
from optparse import OptionParser
import damask

scriptID   = string.replace('$Id$','\n','\\n')
scriptName = os.path.splitext(scriptID.split()[1])[0]

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Average each data block of size 'packing' into single values thus reducing the former grid to grid/packing.

""", version = scriptID)

parser.add_option('-c','--coordinates',
                  dest = 'coords',
                  type = 'string', metavar = 'string',
                  help = 'column heading for coordinates [%default]')
parser.add_option('-p','--packing',
                  dest = 'packing',
                  type = 'int', nargs = 3, metavar = 'int int int',
                  help = 'size of packed group [%default]')
parser.add_option('--shift',
                  dest = 'shift',
                  type = 'int', nargs = 3, metavar = 'int int int',
                  help = 'shift vector of packing stencil [%default]')
parser.add_option('-g', '--grid',
                  dest = 'grid',
                  type = 'int', nargs = 3, metavar = 'int int int',
                  help = 'grid in x,y,z [autodetect]')
parser.add_option('-s', '--size',
                  dest = 'size',
                  type = 'float', nargs = 3, metavar = 'float float float',
                  help = 'size in x,y,z [autodetect]')
parser.set_defaults(coords  = 'ipinitialcoord',
                    packing = (2,2,2),
                    shift   = (0,0,0),
                    grid    = (0,0,0),
                    size    = (0.0,0.0,0.0),
                   )

(options,filenames) = parser.parse_args()

packing = np.array(options.packing,dtype = int)
shift   = np.array(options.shift,  dtype = int)

prefix = 'averagedDown{}x{}x{}_'.format(*packing)
if any(shift != 0): prefix += 'shift{:+}{:+}{:+}_'.format(*shift)

# --- loop over input files ------------------------------------------------------------------------

if filenames == []: filenames = [None]

for name in filenames:
  try:
    table = damask.ASCIItable(name    = name,
                              outname = os.path.join(os.path.dirname(name),
                                                     prefix+os.path.basename(name)) if name else name,
                              buffered = False)
  except: continue
  table.croak('\033[1m'+scriptName+'\033[0m'+(': '+name if name else ''))

# ------------------------------------------ read header ------------------------------------------

  table.head_read()

# ------------------------------------------ sanity checks ----------------------------------------

  errors  = []
  remarks = []
  colCoord = None
  
  if table.label_dimension(options.coords) != 3:  errors.append('coordinates {} are not a vector.'.format(options.coords))
  else: colCoord = table.label_index(options.coords)

  if remarks != []: table.croak(remarks)
  if errors  != []:
    table.croak(errors)
    table.close(dismiss = True)
    continue


# ------------------------------------------ assemble header ---------------------------------------

  table.info_append(scriptID + '\t' + ' '.join(sys.argv[1:]))
  table.head_write()

# --------------- figure out size and grid ---------------------------------------------------------

  table.data_readArray()

  if (any(options.grid) == 0 or any(options.size) == 0.0):
    coords = [np.unique(table.data[:,colCoord+i]) for i in xrange(3)]
    mincorner = np.array(map(min,coords))
    maxcorner = np.array(map(max,coords))
    grid   = np.array(map(len,coords),'i')
    size   = grid/np.maximum(np.ones(3,'d'), grid-1.0) * (maxcorner-mincorner)                        # size from edge to edge = dim * n/(n-1) 
    size   = np.where(grid > 1, size, min(size[grid > 1]/grid[grid > 1]))                             # spacing for grid==1 equal to smallest among other spacings
    delta  = size/np.maximum(np.ones(3,'d'), grid)
    origin = mincorner - 0.5*delta                                                                    # shift from cell center to corner

  else:
    grid   = np.array(options.grid,'i')
    size   = np.array(options.size,'d')
    origin = np.zeros(3,'d')

  packing = np.where(grid == 1,1,packing)                                                           # reset packing to 1 where grid==1
  shift   = np.where(grid == 1,0,shift)                                                             # reset   shift to 0 where grid==1
  packedGrid = np.maximum(np.ones(3,'i'),grid//packing)

  averagedDown = scipy.ndimage.filters.uniform_filter( \
                  np.roll(
                  np.roll(
                  np.roll(table.data.reshape(list(grid)+[table.data.shape[1]],order = 'F'),
                          -shift[0],axis = 0),
                          -shift[1],axis = 1),
                          -shift[2],axis = 2),
                  size = list(packing) + [1],
                  mode = 'wrap',
                  origin = list(-(packing/2)) + [0])\
                  [::packing[0],::packing[1],::packing[2],:].reshape((packedGrid.prod(),table.data.shape[1]),order = 'F')

  
  table.data = averagedDown

#--- generate grid --------------------------------------------------------------------------------

  if colCoord:
    x = (0.5 + shift[0] + np.arange(packedGrid[0],dtype=float))/packedGrid[0]*size[0] + origin[0]
    y = (0.5 + shift[1] + np.arange(packedGrid[1],dtype=float))/packedGrid[1]*size[1] + origin[1]
    z = (0.5 + shift[2] + np.arange(packedGrid[2],dtype=float))/packedGrid[2]*size[2] + origin[2]

    xx = np.tile(          x,                packedGrid[1]* packedGrid[2])
    yy = np.tile(np.repeat(y,packedGrid[0]                ),packedGrid[2])
    zz =         np.repeat(z,packedGrid[0]*packedGrid[1])

    table.data[:,colCoord:colCoord+3] = np.squeeze(np.dstack((xx,yy,zz)))

# ------------------------------------------ output result -----------------------------------------  

  table.data_writeArray()
  
# ------------------------------------------ output finalization -----------------------------------  

  table.close()                                                                                     # close ASCII tables
