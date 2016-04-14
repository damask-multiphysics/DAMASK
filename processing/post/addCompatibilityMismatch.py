#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys
import numpy as np
from optparse import OptionParser
import damask

scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])


def volTetrahedron(coords):
  """
  Return the volume of the tetrahedron with given vertices or sides.
  
  Ifvertices are given they must be in a NumPy array with shape (4,3): the
  position vectors of the 4 vertices in 3 dimensions; if the six sides are
  given, they must be an array of length 6. If both are given, the sides
  will be used in the calculation.

  This method implements
  Tartaglia's formula using the Cayley-Menger determinant:
              |0   1    1    1    1  |
              |1   0   s1^2 s2^2 s3^2|
    288 V^2 = |1  s1^2  0   s4^2 s5^2|
              |1  s2^2 s4^2  0   s6^2|
              |1  s3^2 s5^2 s6^2  0  |
  where s1, s2, ..., s6 are the tetrahedron side lengths.

  from http://codereview.stackexchange.com/questions/77593/calculating-the-volume-of-a-tetrahedron
  """
  # The indexes of rows in the vertices array corresponding to all
  # possible pairs of vertices
  vertex_pair_indexes = np.array(((0, 1), (0, 2), (0, 3),
                                  (1, 2), (1, 3), (2, 3)))

  # Get all the squares of all side lengths from the differences between
  # the 6 different pairs of vertex positions
  vertices =  np.concatenate((coords[0],coords[1],coords[2],coords[3])).reshape([4,3])
  vertex1, vertex2 = vertex_pair_indexes[:,0], vertex_pair_indexes[:,1]
  sides_squared = np.sum((vertices[vertex1] - vertices[vertex2])**2,axis=-1)


  # Set up the Cayley-Menger determinant
  M = np.zeros((5,5))
  # Fill in the upper triangle of the matrix
  M[0,1:] = 1
  # The squared-side length elements can be indexed using the vertex
  # pair indices (compare with the determinant illustrated above)
  M[tuple(zip(*(vertex_pair_indexes + 1)))] = sides_squared

  # The matrix is symmetric, so we can fill in the lower triangle by
  # adding the transpose
  M = M + M.T
  return np.sqrt(np.linalg.det(M) / 288)


def volumeMismatch(size,F,nodes):
  """
  calculates the volume mismatch
  
  volume mismatch is defined as the difference between volume of reconstructed 
  (compatible) cube and determinant of defgrad at the FP
  """
  coords = np.empty([8,3])
  vMismatch = np.empty(grid)
  volInitial = size.prod()/grid.prod()
 
#--------------------------------------------------------------------------------------------------
# calculate actual volume and volume resulting from deformation gradient
  for k in xrange(grid[2]):
    for j in xrange(grid[1]):
      for i in xrange(grid[0]):
        coords[0,0:3] = nodes[0:3,i,  j,  k  ]
        coords[1,0:3] = nodes[0:3,i+1,j,  k  ]
        coords[2,0:3] = nodes[0:3,i+1,j+1,k  ]
        coords[3,0:3] = nodes[0:3,i,  j+1,k  ]
        coords[4,0:3] = nodes[0:3,i,  j,  k+1]
        coords[5,0:3] = nodes[0:3,i+1,j,  k+1]
        coords[6,0:3] = nodes[0:3,i+1,j+1,k+1]
        coords[7,0:3] = nodes[0:3,i,  j+1,k+1]
        vMismatch[i,j,k] = \
           abs(volTetrahedron([coords[6,0:3],coords[0,0:3],coords[7,0:3],coords[3,0:3]])) \
         + abs(volTetrahedron([coords[6,0:3],coords[0,0:3],coords[7,0:3],coords[4,0:3]])) \
         + abs(volTetrahedron([coords[6,0:3],coords[0,0:3],coords[2,0:3],coords[3,0:3]])) \
         + abs(volTetrahedron([coords[6,0:3],coords[0,0:3],coords[2,0:3],coords[1,0:3]])) \
         + abs(volTetrahedron([coords[6,0:3],coords[4,0:3],coords[1,0:3],coords[5,0:3]])) \
         + abs(volTetrahedron([coords[6,0:3],coords[4,0:3],coords[1,0:3],coords[0,0:3]]))
        vMismatch[i,j,k] = vMismatch[i,j,k]/np.linalg.det(F[0:3,0:3,i,j,k])

  return vMismatch/volInitial



def shapeMismatch(size,F,nodes,centres):
  """
  Routine to calculate the shape mismatch
  
  shape mismatch is defined as difference between the vectors from the central point to
  the corners of reconstructed (combatible) volume element and the vectors calculated by deforming
  the initial volume element with the  current deformation gradient
  """
  coordsInitial = np.empty([8,3])
  sMismatch    = np.empty(grid)
   
#--------------------------------------------------------------------------------------------------
# initial positions
  coordsInitial[0,0:3] = [-size[0]/grid[0],-size[1]/grid[1],-size[2]/grid[2]]
  coordsInitial[1,0:3] = [+size[0]/grid[0],-size[1]/grid[1],-size[2]/grid[2]]
  coordsInitial[2,0:3] = [+size[0]/grid[0],+size[1]/grid[1],-size[2]/grid[2]]
  coordsInitial[3,0:3] = [-size[0]/grid[0],+size[1]/grid[1],-size[2]/grid[2]]
  coordsInitial[4,0:3] = [-size[0]/grid[0],-size[1]/grid[1],+size[2]/grid[2]]
  coordsInitial[5,0:3] = [+size[0]/grid[0],-size[1]/grid[1],+size[2]/grid[2]]
  coordsInitial[6,0:3] = [+size[0]/grid[0],+size[1]/grid[1],+size[2]/grid[2]]
  coordsInitial[7,0:3] = [-size[0]/grid[0],+size[1]/grid[1],+size[2]/grid[2]]
  coordsInitial = coordsInitial/2.0
 
#--------------------------------------------------------------------------------------------------
# compare deformed original and deformed positions to actual positions
  for k in xrange(grid[2]):
    for j in xrange(grid[1]):
      for i in xrange(grid[0]):
       sMismatch[i,j,k] = \
         + np.linalg.norm(nodes[0:3,i,  j,    k] - centres[0:3,i,j,k] - np.dot(F[:,:,i,j,k], coordsInitial[0,0:3]))\
         + np.linalg.norm(nodes[0:3,i+1,j,    k] - centres[0:3,i,j,k] - np.dot(F[:,:,i,j,k], coordsInitial[1,0:3]))\
         + np.linalg.norm(nodes[0:3,i+1,j+1,k  ] - centres[0:3,i,j,k] - np.dot(F[:,:,i,j,k], coordsInitial[2,0:3]))\
         + np.linalg.norm(nodes[0:3,i,  j+1,k  ] - centres[0:3,i,j,k] - np.dot(F[:,:,i,j,k], coordsInitial[3,0:3]))\
         + np.linalg.norm(nodes[0:3,i,  j,  k+1] - centres[0:3,i,j,k] - np.dot(F[:,:,i,j,k], coordsInitial[4,0:3]))\
         + np.linalg.norm(nodes[0:3,i+1,j,  k+1] - centres[0:3,i,j,k] - np.dot(F[:,:,i,j,k], coordsInitial[5,0:3]))\
         + np.linalg.norm(nodes[0:3,i+1,j+1,k+1] - centres[0:3,i,j,k] - np.dot(F[:,:,i,j,k], coordsInitial[6,0:3]))\
         + np.linalg.norm(nodes[0:3,i,  j+1,k+1] - centres[0:3,i,j,k] - np.dot(F[:,:,i,j,k], coordsInitial[7,0:3]))
  return sMismatch


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

# ------------------------------------------ assemble header --------------------------------------

  table.info_append(scriptID + '\t' + ' '.join(sys.argv[1:]))
  if options.shape:  table.labels_append('shapeMismatch({})'.format(options.defgrad))
  if options.volume: table.labels_append('volMismatch({})'.format(options.defgrad))

# --------------- figure out size and grid ---------------------------------------------------------

  table.data_readArray()

  coords = [np.unique(table.data[:,colCoord+i]) for i in xrange(3)]
  mincorner = np.array(map(min,coords))
  maxcorner = np.array(map(max,coords))
  grid   = np.array(map(len,coords),'i')
  size   = grid/np.maximum(np.ones(3,'d'), grid-1.0) * (maxcorner-mincorner)                        # size from edge to edge = dim * n/(n-1) 
  size   = np.where(grid > 1, size, min(size[grid > 1]/grid[grid > 1]))                             # spacing for grid==1 set to smallest among other spacings

  N = grid.prod()
  
# --------------- figure out columns to process  ---------------------------------------------------
  key = '1_%s'%options.defgrad
  if key not in table.labels:
    file['croak'].write('column %s not found...\n'%key)
    continue
  else:
    column = table.labels.index(key)                                                                # remember columns of requested data

# ------------------------------------------ assemble header ---------------------------------------
  if options.shape:  table.labels_append(['shapeMismatch(%s)' %options.defgrad])
  if options.volume: table.labels_append(['volMismatch(%s)'%options.defgrad])
  table.head_write()

# ------------------------------------------ read deformation gradient field -----------------------
  table.data_rewind()
  F = np.zeros(N*9,'d').reshape([3,3]+list(grid))
  idx = 0
  while table.data_read():    
    (x,y,z) = damask.util.gridLocation(idx,grid)                                                     # figure out (x,y,z) position from line count
    idx += 1
    F[0:3,0:3,x,y,z] = np.array(map(float,table.data[column:column+9]),'d').reshape(3,3)
  Favg = damask.core.math.tensorAvg(F)
  centres = damask.core.mesh.deformedCoordsFFT(size,F,Favg,[1.0,1.0,1.0])
  
  nodes   = damask.core.mesh.nodesAroundCentres(size,Favg,centres)
  if options.shape:   shapeMismatch = shapeMismatch( size,F,nodes,centres)
  if options.volume: volumeMismatch = volumeMismatch(size,F,nodes)

# ------------------------------------------ process data ------------------------------------------
  table.data_rewind()
  idx = 0
  outputAlive = True
  while outputAlive and table.data_read():                                                          # read next data line of ASCII table
    (x,y,z) = damask.util.gridLocation(idx,grid)                                                    # figure out (x,y,z) position from line count
    idx += 1
    if options.shape:  table.data_append( shapeMismatch[x,y,z])
    if options.volume: table.data_append(volumeMismatch[x,y,z])
    outputAlive = table.data_write()

# ------------------------------------------ output finalization -----------------------------------  

  table.close()                                                                                     # close ASCII tables
