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
Operates on periodic ordered three-dimensional data sets.

""", version = scriptID)

parser.add_option('-c','--coordinates', dest='coords', metavar='string',
                  help='column label of coordinates [%default]')
parser.add_option('-f','--defgrad',     dest='defgrad', metavar='string',
                  help='column label of deformation gradient [%default]')
parser.add_option('--scaling', dest='scaling', type='float', nargs=3, , metavar = ' '.join(['float']*3),
                  help='x/y/z scaling of displacment fluctuation')
parser.set_defaults(coords  = 'ipinitialcoord',
                    defgrad = 'f',
                    scaling = [1.,1.,1.],
                   )

(options,filenames) = parser.parse_args()

# --- loop over input files -------------------------------------------------------------------------
if filenames == []:
  filenames = ['STDIN']

for name in filenames:
  if name == 'STDIN':
    file = {'name':'STDIN', 'input':sys.stdin, 'output':sys.stdout, 'croak':sys.stderr}
    file['croak'].write('\033[1m'+scriptName+'\033[0m\n')
  else:
    if not os.path.exists(name): continue
    file = {'name':name, 'input':open(name), 'output':open(name+'_tmp','w'), 'croak':sys.stderr}
    file['croak'].write('\033[1m'+scriptName+'\033[0m: '+file['name']+'\n')

  table = damask.ASCIItable(file['input'],file['output'],buffered=False)                            # make unbuffered ASCII_table
  table.head_read()                                                                                 # read ASCII header info
  table.info_append(scriptID + '\t' + ' '.join(sys.argv[1:]))

# --------------- figure out columns to process  ---------------------------------------------------

  if table.label_dimension(options.coords) != 3:
    file['croak'].write('no coordinate vector (1/2/3_%s) found...\n'%options.coords)
    continue
  if table.label_dimension(options.defgrad) != 9:
    file['croak'].write('no deformation gradient tensor (1..9_%s) found...\n'%options.defgrad)
    continue

# --------------- figure out size and grid ---------------------------------------------------------

  colCoords  = table.label_index(options.coords)                                                    # starting column of location data
  colDefGrad = table.label_index(options.defgrad)                                                   # remember columns of requested data

  coords = [{},{},{}]
  while table.data_read():                                                                          # read next data line of ASCII table
    for j in xrange(3):
      coords[j][str(table.data[colCoords+j])] = True                                                # remember coordinate along x,y,z
  grid = np.array([len(coords[0]),\
                   len(coords[1]),\
                   len(coords[2]),],'i')                                                            # grid is number of distinct coordinates found
  size = grid/np.maximum(np.ones(3,'d'),grid-1.0)* \
            np.array([max(map(float,coords[0].keys()))-min(map(float,coords[0].keys())),\
                      max(map(float,coords[1].keys()))-min(map(float,coords[1].keys())),\
                      max(map(float,coords[2].keys()))-min(map(float,coords[2].keys())),\
                      ],'d')                                                                        # size from bounding box, corrected for cell-centeredness

  for i, points in enumerate(grid):
    if points == 1:
      options.packing[i] = 1
      options.shift[i]   = 0
      mask = np.ones(3,dtype=bool)
      mask[i]=0
      size[i] = min(size[mask]/grid[mask])                                                          # third spacing equal to smaller of other spacing
  
  N = grid.prod()


# ------------------------------------------ assemble header ---------------------------------------
  table.labels_append(['%s_%s%s'%(coord+1,options.defgrad,options.coords) for coord in xrange(3)])  # extend ASCII header with new labels
  table.head_write()

# ------------------------------------------ read deformation gradient field -----------------------
  table.data_rewind()
  F = np.array([0.0 for i in xrange(N*9)]).reshape([3,3]+list(grid))
  idx = 0
  while table.data_read():
    (x,y,z) = damask.util.gridLocation(idx,grid)                                                    # figure out (x,y,z) position from line count
    idx += 1
    F[0:3,0:3,x,y,z] = np.array(map(float,table.data[colDefGrad:colDefGrad+9]),'d').reshape(3,3)

# ------------------------------------------ calculate coordinates ---------------------------------
  Favg = damask.core.math.tensorAvg(F)
  centroids = damask.core.mesh.deformedCoordsFFT(size,F,Favg,options.scaling)
  
# ------------------------------------------ process data ------------------------------------------
  table.data_rewind()
  idx = 0
  outputAlive = True
  while outputAlive and table.data_read():                                                          # read next data line of ASCII table
    (x,y,z) = damask.util.gridLocation(idx,grid)                                                    # figure out (x,y,z) position from line count
    idx += 1
    table.data_append(list(centroids[:,x,y,z]))
    outputAlive = table.data_write()                                                                # output processed line
  
# ------------------------------------------ output result -----------------------------------------
  outputAlive and table.output_flush()                                                              # just in case of buffered ASCII table

  table.close()                                                                                     # close tables
  os.rename(file['name']+'_tmp',file['name'])                                                       # overwrite old one with tmp new
