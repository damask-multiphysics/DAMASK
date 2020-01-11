#!/usr/bin/env python3

import os
import sys
from io import StringIO
from optparse import OptionParser

import numpy as np
import scipy.ndimage

import damask


scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])


# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [ASCIItable(s)]', description = """
Average each data block of size 'packing' into single values thus reducing the former grid to grid/packing.

""", version = scriptID)

parser.add_option('-c','--coordinates',
                  dest = 'pos',
                  type = 'string', metavar = 'string',
                  help = 'column label of coordinates [%default]')
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
                  help = 'grid in x,y,z (optional)')
parser.add_option('-s', '--size',
                  dest = 'size',
                  type = 'float', nargs = 3, metavar = 'float float float',
                  help = 'size in x,y,z (optional)')
parser.set_defaults(pos     = 'pos',
                    packing = (2,2,2),
                    shift   = (0,0,0),
                   )

(options,filenames) = parser.parse_args()
if filenames == []: filenames = [None]

packing = np.array(options.packing,dtype = int)
shift   = np.array(options.shift,  dtype = int)

prefix = 'averagedDown{}x{}x{}_'.format(*packing)
if any(shift != 0): prefix += 'shift{:+}{:+}{:+}_'.format(*shift)


for name in filenames:
  damask.util.report(scriptName,name)
  
  table = damask.Table.from_ASCII(StringIO(''.join(sys.stdin.read())) if name is None else name)

  if (options.grid is None or options.size is None):
    grid,size,origin = damask.grid_filters.cell_coord0_2_DNA(table.get(options.pos))
  else:
    grid   = np.array(options.grid,'i')
    size   = np.array(options.size,'d')

  packing = np.where(grid == 1,1,packing)                                                           # reset packing to 1 where grid==1
  shift   = np.where(grid == 1,0,shift)                                                             # reset   shift to 0 where grid==1
  packedGrid = np.maximum(np.ones(3,'i'),grid//packing)

  data = table.data.values.reshape(tuple(grid)+(-1,),order = 'F')
  averagedDown = scipy.ndimage.filters.uniform_filter( \
                  np.roll(
                  np.roll(
                  np.roll(data,
                          -shift[0],axis = 0),
                          -shift[1],axis = 1),
                          -shift[2],axis = 2),
                  size = list(packing) + [1],
                  mode = 'wrap',
                  origin = list(-(packing//2)) + [0])\
                  [::packing[0],::packing[1],::packing[2],:].reshape((packedGrid.prod(),-1),order = 'F')

  
  table = damask.Table(averagedDown,table.shapes,table.comments)

  coords = damask.grid_filters.cell_coord0(packedGrid,size,shift/packedGrid*size+origin)
  table.set(options.pos, coords.reshape((-1,3)))


  outname = os.path.join(os.path.dirname(name),prefix+os.path.basename(name))
  table.to_ASCII(sys.stdout if name is None else outname)
