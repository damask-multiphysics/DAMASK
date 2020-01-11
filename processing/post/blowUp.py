#!/usr/bin/env python3

import os
import sys
from io import StringIO
from optparse import OptionParser

from scipy import ndimage
import numpy as np

import damask


scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])


# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [ASCIItable(s)]', description = """
Blows up each value to a surrounding data block of size 'packing' thus increasing the former resolution
to resolution*packing.

""", version = scriptID)

parser.add_option('-c','--coordinates',
                  dest = 'pos', metavar = 'string',
                  help = 'column label of coordinates [%default]')
parser.add_option('-p','--packing',
                  dest = 'packing', type = 'int', nargs = 3, metavar = 'int int int',
                  help = 'dimension of packed group [%default]')
parser.add_option('-g','--grid',
                  dest = 'resolution', type = 'int', nargs = 3, metavar = 'int int int',
                  help = 'grid in x,y,z (optional)')
parser.add_option('-s','--size',
                  dest = 'dimension', type = 'float', nargs = 3, metavar = 'int int int',
                  help = 'size in x,y,z (optional)')
parser.set_defaults(pos  = 'pos',
                    packing = (2,2,2),
                    grid    = (0,0,0),
                    size    = (0.0,0.0,0.0),
                   )

(options,filenames) = parser.parse_args()
if filenames == []: filenames = [None]

options.packing = np.array(options.packing)
prefix = 'blowUp{}x{}x{}_'.format(*options.packing)


for name in filenames:
  damask.util.report(scriptName,name)
  
  table = damask.Table.from_ASCII(StringIO(''.join(sys.stdin.read())) if name is None else name)
  grid,size,origin = damask.grid_filters.cell_coord0_2_DNA(table.get(options.pos))
  
  packing = np.array(options.packing,'i')
  outSize = grid*packing
  
  data = table.data.values.reshape(tuple(grid)+(-1,))
  blownUp = ndimage.interpolation.zoom(data,tuple(packing)+(1,),order=0,mode='nearest').reshape((outSize.prod(),-1))

  table = damask.Table(blownUp,table.shapes,table.comments)

  coords = damask.grid_filters.cell_coord0(outSize,size,origin)
  table.set(options.pos,coords.reshape((-1,3)))
  table.set('elem',np.arange(1,outSize.prod()+1))
  
  outname = os.path.join(os.path.dirname(name),prefix+os.path.basename(name))
  table.to_ASCII(sys.stdout if name is None else outname)
