#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,math,string
import numpy as np
from optparse import OptionParser
import damask

scriptID   = string.replace('$Id: addDeformedConfiguration.py 4500 2015-09-24 09:24:42Z MPIE\m.diehl $','\n','\\n')
scriptName = os.path.splitext(scriptID.split()[1])[0]

#--------------------------------------------------------------------------------------------------
def nodesAroundCentres(gDim,Favg,centres):
#--------------------------------------------------------------------------------------------------
 neighbor = np.array([0, 0, 0,
                      1, 0, 0,
                      1, 1, 0,
                      0, 1, 0,
                      0, 0, 1,
                      1, 0, 1,
                      1, 1, 1,
                      0, 1, 1]).reshape(8,3)

#--------------------------------------------------------------------------------------------------
# building wrappedCentres = centroids + ghosts
 diag = np.ones([3])
 wrappedCentres = np.zeros([3,grid[0]+2,grid[1]+2,grid[2]+2])
 wrappedCentres[0:3,1:grid[0]+1,1:grid[1]+1,1:grid[2]+1] = centres
 for k in xrange(grid[2]+2):
   for j in xrange(grid[1]+2):
     for i in xrange(grid[0]+2):
       if (k in [0,grid[2]+1] or j in [0,grid[1]+1] or i in[0,grid[0]+1]):
         me = np.array([i,j,k],'i')                                                                  # me on skin
         shift = abs(grid+np.ones([3],'i')-2*me)/(grid+np.ones([3],'i'))*\
                                                                np.sign(grid+np.ones([3],'i')-2*me)
         lookup = np.array(me-diag+shift*grid,'i')
         wrappedCentres[0:3,i,        j,        k] = \
                centres[0:3,lookup[0],lookup[1],lookup[2]] - np.dot(Favg, shift*gDim)
 
#--------------------------------------------------------------------------------------------------
# averaging
 nodes = np.zeros([3,grid[0]+1,grid[1]+1,grid[2]+1])
 for k in xrange(grid[2]+1):
   for j in xrange(grid[1]+1):
     for i in xrange(grid[0]+1):
       for n in xrange(8):
         nodes[0:3,i,j,k] = \
         nodes[0:3,i,j,k] + wrappedCentres[0:3,i+neighbor[n,0],j+neighbor[n,1],k+neighbor[n,2] ]

 return nodes/8.0

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options file[s]', description = """
Add deformed configuration of given initial coordinates.
Operates on periodic three-dimensional x,y,z-ordered data sets.

""", version = scriptID)

parser.add_option('-f', '--defgrad',dest='defgrad', metavar = 'string',
                                    help='heading of deformation gradient columns [%default]')
parser.add_option('-u', '--unitlength', dest='unitlength', type='float', metavar = 'float',
                                    help='set unit length for 2D model [%default]')

parser.set_defaults(deformed = 'ipinitialcoord')
parser.set_defaults(unitlength = 0.0)

(options,filenames) = parser.parse_args()

options.scaling += [1.0 for i in xrange(max(0,3-len(options.scaling)))]
scaling = map(float, options.scaling)


# --- loop over input files -------------------------------------------------------------------------

if filenames == []: filenames = [None]

for name in filenames:
  try:
    table = damask.ASCIItable(name = name,
                              buffered = False)
  except: continue
  damask.util.report(scriptName,name)

# ------------------------------------------ read header ------------------------------------------

  table.head_read()

# ------------------------------------------ sanity checks ----------------------------------------

  errors  = []
  remarks = []
  
  if table.label_dimension(options.coords) != 3:  errors.append('coordinates {} are not a vector.'.format(options.coords))
  else: colCoord = table.label_index(options.coords)

  if table.label_dimension(options.defgrad) != 9: errors.append('deformation gradient {} is not a tensor.'.format(options.defgrad))
  else: colF = table.label_index(options.defgrad)

  if remarks != []: damask.util.croak(remarks)
  if errors  != []:
    damask.util.croak(errors)
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
    damask.util.croak(errors)
    table.close(dismiss = True)
    continue
  
# ------------------------------------------ assemble header ---------------------------------------

  table.info_append(scriptID + '\t' + ' '.join(sys.argv[1:]))
  for coord in xrange(3):
    label = '{}_{}_{}'.format(coord+1,options.defgrad,options.coords)
    if np.any(scaling) != 1.0: label+='_{}_{}_{}'.format(scaling)
    if options.undeformed: label+='_undeformed'
    table.labels_append([label])                                                                    # extend ASCII header with new labels
  table.head_write()

# ------------------------------------------ read deformation gradient field -----------------------
  centroids,Favg = deformedCoordsFFT(table.data[:,colF:colF+9].reshape(grid[0],grid[1],grid[2],3,3))

# ------------------------------------------ process data ------------------------------------------
  table.data_rewind()
  idx = 0
  outputAlive = True
  while outputAlive and table.data_read():                                                          # read next data line of ASCII table
    (x,y,z) = damask.util.gridLocation(idx,grid)                                                    # figure out (x,y,z) position from line count
    idx += 1
    table.data_append(list(centroids[z,y,x,:]))
    outputAlive = table.data_write()                       

# ------------------------------------------ output finalization -----------------------------------  

  table.close()                                                                                     # close ASCII tables
