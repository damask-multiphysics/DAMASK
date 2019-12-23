#!/usr/bin/env python3

import os
import sys
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

options.packing = np.array(options.packing)
prefix = 'blowUp{}x{}x{}_'.format(*options.packing)

# --- loop over input files -------------------------------------------------------------------------

if filenames == []: filenames = [None]

for name in filenames:
  try:    table = damask.ASCIItable(name = name,
                                    outname = os.path.join(os.path.dirname(name),
                                                           prefix+os.path.basename(name)) if name else name,
                                    buffered = False)
  except IOError:
    continue
  damask.util.report(scriptName,name)

# ------------------------------------------ read header ------------------------------------------

  table.head_read()

# ------------------------------------------ sanity checks ----------------------------------------
  
  errors  = []
  remarks = []
  
  if table.label_dimension(options.pos) != 3:  errors.append('coordinates "{}" are not a vector.'.format(options.pos))

  colElem = table.label_index('elem')
  
  if remarks != []: damask.util.croak(remarks)
  if errors  != []:
    damask.util.croak(errors)
    table.close(dismiss = True)
    continue

# --------------- figure out size and grid ---------------------------------------------------------

  table.data_readArray(options.pos)
  table.data_rewind()

  grid,size,origin = damask.grid_filters.cell_coord0_2_DNA(table.data)
  
  packing = np.array(options.packing,'i')
  outSize = grid*packing
  
# ------------------------------------------ assemble header --------------------------------------

  table.info_append(scriptID + '\t' + ' '.join(sys.argv[1:]))
  table.head_write()

# ------------------------------------------ process data -------------------------------------------
  table.data_readArray()
  data = table.data.reshape(tuple(grid)+(-1,))
  table.data = ndimage.interpolation.zoom(data,tuple(packing)+(1,),order=0,mode='nearest').reshape((outSize.prod(),-1))
  coords = damask.grid_filters.cell_coord0(outSize,size,origin)
  table.data[:,table.label_indexrange(options.pos)] = coords.reshape((-1,3))
  table.data[:,table.label_index('elem')] = np.arange(1,outSize.prod()+1)

# ------------------------------------------ output finalization -----------------------------------  
  table.data_writeArray()
  table.close()
