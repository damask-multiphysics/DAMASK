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

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Blows up each value to a surrounding data block of size 'packing' thus increasing the former resolution
to resolution*packing.

""", version = scriptID)

parser.add_option('-c','--coordinates', dest='coords', metavar='string',
                  help='column heading for coordinates [%default]')
parser.add_option('-p','--packing',     dest='packing', type='int', nargs=3, metavar='int int int',
                  help='dimension of packed group [%default]')
parser.add_option('-g','--grid',        dest='resolution', type='int', nargs=3, metavar='int int int',
                  help='resolution in x,y,z [autodetect]')
parser.add_option('-s','--size',        dest='dimension', type='float', nargs=3, metavar='int int int',
                  help='dimension in x,y,z [autodetect]')
parser.set_defaults(coords  = 'ipinitialcoord')
parser.set_defaults(packing = (2,2,2))
parser.set_defaults(grid    = (0,0,0))
parser.set_defaults(size    = (0.0,0.0,0.0))

(options,filenames) = parser.parse_args()


options.packing = np.array(options.packing)
prefix = 'blowUp%ix%ix%i_'%(options.packing[0],options.packing[1],options.packing[2])

# --- loop over input files -------------------------------------------------------------------------

if filenames == []: filenames = [None]

for name in filenames:
  try:
    table = damask.ASCIItable(name = name,
                              outname = os.path.join(os.path.dirname(name),
                                                     prefix+ \
                                                     os.path.basename(name)) if name else name,
                              buffered = False)
  except: continue
  table.croak(damask.util.emph(scriptName)+(': '+name if name else ''))

# ------------------------------------------ read header ------------------------------------------

  table.head_read()

# ------------------------------------------ sanity checks ----------------------------------------
  
  if table.label_dimension(options.coords) != 3: errors.append('coordinates {} are not a vector.'.format(options.coords))
  else: coordCol = table.label_index(options.coords)


# ------------------------------------------ assemble header --------------------------------------

  table.info_append(scriptID + '\t' + ' '.join(sys.argv[1:]))

# --------------- figure out size and grid ---------------------------------------------------------

  table.data_readArray()

  coords = [{},{},{}]
  for i in xrange(len(table.data)):  
    for j in xrange(3):
      coords[j][str(table.data[i,coordCol+j])] = True
  grid = np.array(map(len,coords),'i')
  size = grid/np.maximum(np.ones(3,'d'),grid-1.0)* \
            np.array([max(map(float,coords[0].keys()))-min(map(float,coords[0].keys())),\
                      max(map(float,coords[1].keys()))-min(map(float,coords[1].keys())),\
                      max(map(float,coords[2].keys()))-min(map(float,coords[2].keys())),\
                      ],'d')                                                                        # size from bounding box, corrected for cell-centeredness

  size = np.where(grid > 1, size, min(size[grid > 1]/grid[grid > 1]))                               # spacing for grid==1 equal to smallest among other spacings

  
  packing = np.array(options.packing,'i')
  outSize = grid*packing
  
# ------------------------------------------ assemble header ---------------------------------------
  table.head_write()

# ------------------------------------------ process data -------------------------------------------
  table.data_rewind()
  data = np.zeros(outSize.tolist()+[len(table.labels)])
  p = np.zeros(3,'i')
  
  for p[2] in xrange(grid[2]):
    for p[1] in xrange(grid[1]):
      for p[0] in xrange(grid[0]):
        d = p*packing
        table.data_read()
        data[d[0]:d[0]+packing[0],
             d[1]:d[1]+packing[1],
              d[2]:d[2]+packing[2],
             : ] = np.tile(np.array(table.data_asFloat(),'d'),packing.tolist()+[1])                 # tile to match blowUp voxel size
  elementSize = size/grid/packing
  elem = 1
  for c in xrange(outSize[2]):
    for b in xrange(outSize[1]):
      for a in xrange(outSize[0]):
        data[a,b,c,coordCol:coordCol+3] = [a+0.5,b+0.5,c+0.5]*elementSize
        data[a,b,c,table.label_index('elem')] = elem
        table.data = data[a,b,c,:].tolist()
        outputAlive = table.data_write()                                                            # output processed line
        elem += 1

# ------------------------------------------ output finalization -----------------------------------

  table.close()                                                                                     # close input ASCII table (works for stdin)
