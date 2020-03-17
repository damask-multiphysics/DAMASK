#!/usr/bin/env python3

import os
import sys
from io import StringIO
from optparse import OptionParser

import numpy as np

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
  vertices =  np.concatenate((coords[0],coords[1],coords[2],coords[3])).reshape(4,3)
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
  Calculates the volume mismatch.
  
  volume mismatch is defined as the difference between volume of reconstructed 
  (compatible) cube and determinant of deformation gradient at Fourier point.
  """
  coords = np.empty([8,3])
  vMismatch = np.empty(grid[::-1])
  volInitial = size.prod()/grid.prod()
 
#--------------------------------------------------------------------------------------------------
# calculate actual volume and volume resulting from deformation gradient
  for k in range(grid[2]):
    for j in range(grid[1]):
      for i in range(grid[0]):
        coords[0,0:3] = nodes[k,  j,  i  ,0:3]
        coords[1,0:3] = nodes[k  ,j,  i+1,0:3]
        coords[2,0:3] = nodes[k  ,j+1,i+1,0:3]
        coords[3,0:3] = nodes[k,  j+1,i  ,0:3]
        coords[4,0:3] = nodes[k+1,j,  i  ,0:3]
        coords[5,0:3] = nodes[k+1,j,  i+1,0:3]
        coords[6,0:3] = nodes[k+1,j+1,i+1,0:3]
        coords[7,0:3] = nodes[k+1,j+1,i  ,0:3]
        vMismatch[k,j,i] = \
        (  abs(volTetrahedron([coords[6,0:3],coords[0,0:3],coords[7,0:3],coords[3,0:3]])) \
         + abs(volTetrahedron([coords[6,0:3],coords[0,0:3],coords[7,0:3],coords[4,0:3]])) \
         + abs(volTetrahedron([coords[6,0:3],coords[0,0:3],coords[2,0:3],coords[3,0:3]])) \
         + abs(volTetrahedron([coords[6,0:3],coords[0,0:3],coords[2,0:3],coords[1,0:3]])) \
         + abs(volTetrahedron([coords[6,0:3],coords[4,0:3],coords[1,0:3],coords[5,0:3]])) \
         + abs(volTetrahedron([coords[6,0:3],coords[4,0:3],coords[1,0:3],coords[0,0:3]]))) \
        /np.linalg.det(F[k,j,i,0:3,0:3])
  return vMismatch/volInitial



def shapeMismatch(size,F,nodes,centres):
  """
  Routine to calculate the shape mismatch.
  
  shape mismatch is defined as difference between the vectors from the central point to
  the corners of reconstructed (combatible) volume element and the vectors calculated by deforming
  the initial volume element with the  current deformation gradient.
  """
  coordsInitial = np.empty([8,3])
  sMismatch    = np.empty(grid[::-1])
   
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
  for k in range(grid[2]):
    for j in range(grid[1]):
      for i in range(grid[0]):
       sMismatch[k,j,i] = \
         + np.linalg.norm(nodes[k,  j,  i  ,0:3] - centres[k,j,i,0:3] - np.dot(F[k,j,i,:,:], coordsInitial[0,0:3]))\
         + np.linalg.norm(nodes[k,  j,  i+1,0:3] - centres[k,j,i,0:3] - np.dot(F[k,j,i,:,:], coordsInitial[1,0:3]))\
         + np.linalg.norm(nodes[k,  j+1,i+1,0:3] - centres[k,j,i,0:3] - np.dot(F[k,j,i,:,:], coordsInitial[2,0:3]))\
         + np.linalg.norm(nodes[k,  j+1,i  ,0:3] - centres[k,j,i,0:3] - np.dot(F[k,j,i,:,:], coordsInitial[3,0:3]))\
         + np.linalg.norm(nodes[k+1,j,  i  ,0:3] - centres[k,j,i,0:3] - np.dot(F[k,j,i,:,:], coordsInitial[4,0:3]))\
         + np.linalg.norm(nodes[k+1,j,  i+1,0:3] - centres[k,j,i,0:3] - np.dot(F[k,j,i,:,:], coordsInitial[5,0:3]))\
         + np.linalg.norm(nodes[k+1,j+1,i+1,0:3] - centres[k,j,i,0:3] - np.dot(F[k,j,i,:,:], coordsInitial[6,0:3]))\
         + np.linalg.norm(nodes[k+1,j+1,i  ,0:3] - centres[k,j,i,0:3] - np.dot(F[k,j,i,:,:], coordsInitial[7,0:3]))
  return sMismatch


# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [ASCIItable(s)]', description = """
Add column(s) containing the shape and volume mismatch resulting from given deformation gradient.
Operates on periodic three-dimensional x,y,z-ordered data sets.

""", version = scriptID)


parser.add_option('-c','--coordinates',
                  dest = 'pos',
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
parser.set_defaults(pos     = 'pos',
                    defgrad = 'f',
                    shape   = True,
                    volume  = True,
                   )

(options,filenames) = parser.parse_args()
if filenames == []: filenames = [None]


for name in filenames:
  damask.util.report(scriptName,name)
  
  table = damask.Table.from_ASCII(StringIO(''.join(sys.stdin.read())) if name is None else name)
  grid,size,origin = damask.grid_filters.cell_coord0_gridSizeOrigin(table.get(options.pos))

  F = table.get(options.defgrad).reshape(grid[2],grid[1],grid[0],3,3)
  nodes = damask.grid_filters.node_coord(size,F)
  
  if options.shape:
    centers = damask.grid_filters.cell_coord(size,F)
    shapeMismatch = shapeMismatch(  size,table.get(options.defgrad).reshape(grid[2],grid[1],grid[0],3,3),nodes,centers)
    table.add('shapeMismatch(({}))'.format(options.defgrad),
              shapeMismatch.reshape(-1,1),
              scriptID+' '+' '.join(sys.argv[1:]))
    
  if options.volume:
    volumeMismatch = volumeMismatch(size,table.get(options.defgrad).reshape(grid[2],grid[1],grid[0],3,3),nodes)
    table.add('volMismatch(({}))'.format(options.defgrad),
              volumeMismatch.reshape(-1,1),
              scriptID+' '+' '.join(sys.argv[1:]))

  table.to_ASCII(sys.stdout if name is None else name)
