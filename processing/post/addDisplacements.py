#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,math
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
  
  for i in xrange(datalen):
    node = scipy.ndimage.convolve(cellData.reshape(tuple(grid)+(datalen,))[...,i],
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
def displacementAvgFFT(F,grid,size,nodal=False,transformed=False):
  """calculate average cell center (or nodal) displacement for deformation gradient field specified in each grid cell"""
  if nodal:
    x, y, z = np.meshgrid(np.linspace(0,size[0],1+grid[0]),
                          np.linspace(0,size[1],1+grid[1]),
                          np.linspace(0,size[2],1+grid[2]),
                          indexing = 'ij')
  else:
    x, y, z = np.meshgrid(np.linspace(0,size[0],grid[0],endpoint=False),
                          np.linspace(0,size[1],grid[1],endpoint=False),
                          np.linspace(0,size[2],grid[2],endpoint=False),
                          indexing = 'ij')

  origCoords = np.concatenate((z[:,:,:,None],y[:,:,:,None],x[:,:,:,None]),axis = 3) 

  F_fourier = F if transformed else np.fft.rfftn(F,axes=(0,1,2))                                    # transform or use provided data
  Favg = np.real(F_fourier[0,0,0,:,:])/grid.prod()                                                  # take zero freq for average
  avgDisplacement = np.einsum('ml,ijkl->ijkm',Favg-np.eye(3),origCoords)                            # dX = Favg.X

  return avgDisplacement

#--------------------------------------------------------------------------------------------------
def displacementFluctFFT(F,grid,size,nodal=False,transformed=False):
  """calculate cell center (or nodal) displacement for deformation gradient field specified in each grid cell"""
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

  displacement = np.fft.irfftn(displacement_fourier,grid,axes=(0,1,2))

  return cell2node(displacement,grid) if nodal else displacement


# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options file[s]', description = """
Add displacments resulting from deformation gradient field.
Operates on periodic three-dimensional x,y,z-ordered data sets.
Outputs at cell centers or cell nodes (into separate file).

""", version = scriptID)

parser.add_option('-f', '--defgrad',
                  dest    = 'defgrad',
                  metavar = 'string',
                  help    = 'column label of deformation gradient [%default]')
parser.add_option('-c', '--coordinates',
                  dest    = 'coords',
                  metavar = 'string',
                  help    = 'column label of coordinates [%default]')
parser.add_option('--nodal',
                  dest    = 'nodal',
                  action  = 'store_true',
                  help    = 'output nodal (not cell-centered) displacements')

parser.set_defaults(defgrad = 'f',
                    coords  = 'pos',
                    nodal   = False,
                   )

(options,filenames) = parser.parse_args()

# --- loop over input files -------------------------------------------------------------------------

if filenames == []: filenames = [None]

for name in filenames:
  try:    table = damask.ASCIItable(name = name,
                                    outname = (os.path.splitext(name)[0]+
                                               '_nodal'+
                                               os.path.splitext(name)[1]) if (options.nodal and name) else None,
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

  coordDim = table.label_dimension(options.coords)
  if not 3 >= coordDim >= 1:
    errors.append('coordinates "{}" need to have one, two, or three dimensions.'.format(options.coords))
  elif coordDim < 3:
    remarks.append('appending {} dimension{} to coordinates "{}"...'.format(3-coordDim,
                                                                            's' if coordDim < 2 else '',
                                                                            options.coords))

  if remarks != []: damask.util.croak(remarks)
  if errors  != []:
    damask.util.croak(errors)
    table.close(dismiss=True)
    continue

# --------------- figure out size and grid ---------------------------------------------------------

  table.data_readArray([options.defgrad,options.coords])
  table.data_rewind()

  if len(table.data.shape) < 2: table.data.shape += (1,)                                            # expand to 2D shape
  if table.data[:,9:].shape[1] < 3:
    table.data = np.hstack((table.data,
                            np.zeros((table.data.shape[0],
                                      3-table.data[:,9:].shape[1]),dtype='f')))                     # fill coords up to 3D with zeros

  if remarks != []: damask.util.croak(remarks)
  if errors  != []:
    damask.util.croak(errors)
    table.close(dismiss = True)
    continue

# --------------- figure out size and grid ---------------------------------------------------------

  coords = [np.unique(table.data[:,9+i]) for i in xrange(3)]
  mincorner = np.array(map(min,coords))
  maxcorner = np.array(map(max,coords))
  grid   = np.array(map(len,coords),'i')
  size   = grid/np.maximum(np.ones(3,'d'), grid-1.0) * (maxcorner-mincorner)                        # size from edge to edge = dim * n/(n-1) 
  size   = np.where(grid > 1, size, min(size[grid > 1]/grid[grid > 1]))                             # spacing for grid==1 set to smallest among other spacings

  N = grid.prod()

  if N != len(table.data): errors.append('data count {} does not match grid {}x{}x{}.'.format(N,*grid))
  if errors  != []:
    damask.util.croak(errors)
    table.close(dismiss = True)
    continue
  
# ------------------------------------------ process data ------------------------------------------

  F_fourier = np.fft.rfftn(table.data[:,:9].reshape(grid[2],grid[1],grid[0],3,3),axes=(0,1,2))      # perform transform only once...

  displacement    = displacementFluctFFT(F_fourier,grid,size,options.nodal,transformed=True)
  avgDisplacement = displacementAvgFFT  (F_fourier,grid,size,options.nodal,transformed=True)

# ------------------------------------------ assemble header ---------------------------------------

  if options.nodal:
    table.info_clear()
    table.labels_clear()

  table.info_append(scriptID + '\t' + ' '.join(sys.argv[1:]))
  table.labels_append((['{}_pos'         .format(i+1)      for i in xrange(3)] if options.nodal else []) +
                       ['{}_avg({}).{}'  .format(i+1,options.defgrad,options.coords) for i in xrange(3)] +
                       ['{}_fluct({}).{}'.format(i+1,options.defgrad,options.coords) for i in xrange(3)] )
  table.head_write()

# ------------------------------------------ output data -------------------------------------------

  zrange = np.linspace(0,size[2],1+grid[2]) if options.nodal else xrange(grid[2])
  yrange = np.linspace(0,size[1],1+grid[1]) if options.nodal else xrange(grid[1])
  xrange = np.linspace(0,size[0],1+grid[0]) if options.nodal else xrange(grid[0])

  for i,z     in enumerate(zrange):
    for j,y   in enumerate(yrange):
      for k,x in enumerate(xrange):
        if options.nodal: table.data_clear()
        else:             table.data_read()
        table.data_append([x,y,z] if options.nodal else [])
        table.data_append(list(avgDisplacement[i,j,k,:]))
        table.data_append(list(   displacement[i,j,k,:]))
        table.data_write()                       

# ------------------------------------------ output finalization -----------------------------------

  table.close()                                                                                     # close ASCII tables
