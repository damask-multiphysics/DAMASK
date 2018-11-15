#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os
import math
import numpy as np
import scipy.ndimage
from optparse import OptionParser
import damask

scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])

#--------------------------------------------------------------------------------------------------
def cell2node(cellData,grid):

  nodeData = 0.0
  datalen = np.array(cellData.shape[3:]).prod()
  
  for i in range(datalen):
    node = scipy.ndimage.convolve(cellData.reshape(tuple(grid[::-1])+(datalen,))[...,i],
                                  np.ones((2,2,2))/8.,                                              # 2x2x2 neighborhood of cells
                                  mode = 'wrap',
                                  origin = -1,                                                      # offset to have cell origin as center
                                 )                                                                  # now averaged at cell origins
    node = np.append(node,node[np.newaxis,0,:,:,...],axis=0)                                        # wrap along z
    node = np.append(node,node[:,0,np.newaxis,:,...],axis=1)                                        # wrap along y
    node = np.append(node,node[:,:,0,np.newaxis,...],axis=2)                                        # wrap along x

    nodeData = node[...,np.newaxis] if i==0 else np.concatenate((nodeData,node[...,np.newaxis]),axis=-1)

  return nodeData

#--------------------------------------------------------------------------------------------------
def deformationAvgFFT(F,grid,size,nodal=False,transformed=False):
  """Calculate average cell center (or nodal) deformation for deformation gradient field specified in each grid cell"""
  if nodal:
    x, y, z = np.meshgrid(np.linspace(0,size[2],1+grid[2]),
                          np.linspace(0,size[1],1+grid[1]),
                          np.linspace(0,size[0],1+grid[0]),
                          indexing = 'ij')
  else:
    x, y, z = np.meshgrid(np.linspace(size[2]/grid[2]/2.,size[2]-size[2]/grid[2]/2.,grid[2]),
                          np.linspace(size[1]/grid[1]/2.,size[1]-size[1]/grid[1]/2.,grid[1]),
                          np.linspace(size[0]/grid[0]/2.,size[0]-size[0]/grid[0]/2.,grid[0]),
                          indexing = 'ij')

  origCoords = np.concatenate((z[:,:,:,None],y[:,:,:,None],x[:,:,:,None]),axis = 3) 

  F_fourier = F if transformed else np.fft.rfftn(F,axes=(0,1,2))                                    # transform or use provided data
  Favg = np.real(F_fourier[0,0,0,:,:])/grid.prod()                                                  # take zero freq for average
  avgDeformation = np.einsum('ml,ijkl->ijkm',Favg,origCoords)                                       # dX = Favg.X

  return avgDeformation

#--------------------------------------------------------------------------------------------------
def displacementFluctFFT(F,grid,size,nodal=False,transformed=False):
  """Calculate cell center (or nodal) displacement for deformation gradient field specified in each grid cell"""
  integrator = 0.5j * size / math.pi

  kk, kj, ki = np.meshgrid(np.where(np.arange(grid[2])>grid[2]//2,np.arange(grid[2])-grid[2],np.arange(grid[2])),
                           np.where(np.arange(grid[1])>grid[1]//2,np.arange(grid[1])-grid[1],np.arange(grid[1])),
                                    np.arange(grid[0]//2+1),
                           indexing = 'ij')
  k_s = np.concatenate((ki[:,:,:,None],kj[:,:,:,None],kk[:,:,:,None]),axis = 3) 
  k_sSquared = np.einsum('...l,...l',k_s,k_s)
  k_sSquared[0,0,0] = 1.0                                                                           # ignore global average frequency

#--------------------------------------------------------------------------------------------------
# integration in Fourier space

  displacement_fourier = -np.einsum('ijkml,ijkl,l->ijkm',
                                    F if transformed else np.fft.rfftn(F,axes=(0,1,2)),
                                    k_s,
                                    integrator,
                                    ) / k_sSquared[...,np.newaxis]

#--------------------------------------------------------------------------------------------------
# backtransformation to real space

  displacement = np.fft.irfftn(displacement_fourier,grid[::-1],axes=(0,1,2))

  return cell2node(displacement,grid) if nodal else displacement


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
  Calculates the volume mismatch
  
  volume mismatch is defined as the difference between volume of reconstructed 
  (compatible) cube and determinant of defgrad at the FP
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
  Routine to calculate the shape mismatch
  
  shape mismatch is defined as difference between the vectors from the central point to
  the corners of reconstructed (combatible) volume element and the vectors calculated by deforming
  the initial volume element with the  current deformation gradient
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

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options file[s]', description = """
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
  
  if table.label_dimension(options.defgrad) != 9:
    errors.append('deformation gradient "{}" is not a 3x3 tensor.'.format(options.defgrad))

  coordDim = table.label_dimension(options.pos)
  if not 3 >= coordDim >= 1:
    errors.append('coordinates "{}" need to have one, two, or three dimensions.'.format(options.pos))
  elif coordDim < 3:
    remarks.append('appending {} dimension{} to coordinates "{}"...'.format(3-coordDim,
                                                                            's' if coordDim < 2 else '',
                                                                            options.pos))

  if remarks != []: damask.util.croak(remarks)
  if errors  != []:
    damask.util.croak(errors)
    table.close(dismiss=True)
    continue

# --------------- figure out size and grid ---------------------------------------------------------

  table.data_readArray([options.defgrad,options.pos])
  table.data_rewind()

  if table.data[:,9:].shape[1] < 3:
    table.data = np.hstack((table.data,
                            np.zeros((table.data.shape[0],
                                      3-table.data[:,9:].shape[1]),dtype='f')))                     # fill coords up to 3D with zeros

  grid,size = damask.util.coordGridAndSize(table.data[:,9:12])
  N = grid.prod()

  if N != len(table.data): errors.append('data count {} does not match grid {}x{}x{}.'.format(N,*grid))
  if errors  != []:
    damask.util.croak(errors)
    table.close(dismiss = True)
    continue
  
# -----------------------------process data and assemble header -------------------------------------

  F_fourier = np.fft.rfftn(table.data[:,:9].reshape(grid[2],grid[1],grid[0],3,3),axes=(0,1,2))      # perform transform only once...
  nodes = displacementFluctFFT(F_fourier,grid,size,True,transformed=True)\
        + deformationAvgFFT   (F_fourier,grid,size,True,transformed=True)

  if options.shape:
    table.labels_append(['shapeMismatch({})'.format(options.defgrad)])
    centres = displacementFluctFFT(F_fourier,grid,size,False,transformed=True)\
            + deformationAvgFFT   (F_fourier,grid,size,False,transformed=True)

  if options.volume:
    table.labels_append(['volMismatch({})'.format(options.defgrad)])

  table.head_write()
  if options.shape:
    shapeMismatch = shapeMismatch( size,table.data[:,:9].reshape(grid[2],grid[1],grid[0],3,3),nodes,centres)
  if options.volume:
    volumeMismatch = volumeMismatch(size,table.data[:,:9].reshape(grid[2],grid[1],grid[0],3,3),nodes)

# ------------------------------------------ output data -------------------------------------------
  for i in range(grid[2]):
    for j in range(grid[1]):
      for k in range(grid[0]):
        table.data_read()
        if options.shape:  table.data_append(shapeMismatch[i,j,k])
        if options.volume: table.data_append(volumeMismatch[i,j,k])
        table.data_write()

# ------------------------------------------ output finalization -----------------------------------

  table.close()                                                                                     # close ASCII tables
