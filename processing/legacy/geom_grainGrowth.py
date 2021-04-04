#!/usr/bin/env python3

import os
import sys
from io import StringIO
from optparse import OptionParser

import numpy as np
from scipy import ndimage

import damask


scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])


getInterfaceEnergy = lambda A,B: np.float32((A != B)*1.0)                                           # 1.0 if A & B are distinct, 0.0 otherwise
struc = ndimage.generate_binary_structure(3,1)                                                      # 3D von Neumann neighborhood


#--------------------------------------------------------------------------------------------------
#                                MAIN
#--------------------------------------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog option(s) [geomfile(s)]', description = """
Smoothen interface roughness by simulated curvature flow.
This is achieved by the diffusion of each initially sharply bounded grain volume within the periodic domain
up to a given distance 'd' voxels.
The final geometry is assembled by selecting at each voxel that grain index for which the concentration remains largest.

""", version = scriptID)

parser.add_option('-d', '--distance',
                  dest = 'd',
                  type = 'float', metavar = 'float',
                  help = 'diffusion distance in voxels [%default]')
parser.add_option('-N', '--iterations',
                  dest = 'N',
                  type = 'int', metavar = 'int',
                  help = 'curvature flow iterations [%default]')
parser.add_option('-i', '--immutable',
                  action = 'extend', dest = 'immutable', metavar = '<int LIST>',
                  help = 'list of immutable material indices')
parser.add_option('--ndimage',
                  dest = 'ndimage', action='store_true',
                  help = 'use ndimage.gaussian_filter in lieu of explicit FFT')

parser.set_defaults(d = 1,
                    N = 1,
                    immutable = [],
                    ndimage = False,
                   )

(options, filenames) = parser.parse_args()

options.immutable = list(map(int,options.immutable))


if filenames == []: filenames = [None]

for name in filenames:
  damask.util.report(scriptName,name)

  geom = damask.Grid.load(StringIO(''.join(sys.stdin.read())) if name is None else name)

  grid_original = geom.cells
  damask.util.croak(geom)
  material = np.tile(geom.material,np.where(grid_original == 1, 2,1))                   # make one copy along dimensions with grid == 1
  grid = np.array(material.shape)

# --- initialize support data ---------------------------------------------------------------------

# store a copy of the initial material indices to find locations of immutable indices
  material_original = np.copy(material)

  if not options.ndimage:
    X,Y,Z = np.mgrid[0:grid[0],0:grid[1],0:grid[2]]

    # Calculates gaussian weights for simulating 3d diffusion
    gauss = np.exp(-(X*X + Y*Y + Z*Z)/(2.0*options.d*options.d),dtype=np.float32) \
            /np.power(2.0*np.pi*options.d*options.d,(3.0 - np.count_nonzero(grid_original == 1))/2.,dtype=np.float32)

    gauss[:,:,:grid[2]//2:-1] = gauss[:,:,1:(grid[2]+1)//2]     # trying to cope with uneven (odd) grid size
    gauss[:,:grid[1]//2:-1,:] = gauss[:,1:(grid[1]+1)//2,:]
    gauss[:grid[0]//2:-1,:,:] = gauss[1:(grid[0]+1)//2,:,:]
    gauss = np.fft.rfftn(gauss).astype(np.complex64)

  for smoothIter in range(options.N):

    interfaceEnergy = np.zeros(material.shape,dtype=np.float32)
    for i in (-1,0,1):
      for j in (-1,0,1):
        for k in (-1,0,1):
          # assign interfacial energy to all voxels that have a differing neighbor (in Moore neighborhood)
          interfaceEnergy = np.maximum(interfaceEnergy,
                                       getInterfaceEnergy(material,np.roll(np.roll(np.roll(
                                                          material,i,axis=0), j,axis=1), k,axis=2)))

    # periodically extend interfacial energy array by half a grid size in positive and negative directions
    periodic_interfaceEnergy = np.tile(interfaceEnergy,(3,3,3))[grid[0]//2:-grid[0]//2,
                                                                grid[1]//2:-grid[1]//2,
                                                                grid[2]//2:-grid[2]//2]

    # transform bulk volume (i.e. where interfacial energy remained zero), store index of closest boundary voxel
    index = ndimage.morphology.distance_transform_edt(periodic_interfaceEnergy == 0.,
                                                      return_distances = False,
                                                      return_indices = True)

    # want array index of nearest voxel on periodically extended boundary
    periodic_bulkEnergy = periodic_interfaceEnergy[index[0],
                                                   index[1],
                                                   index[2]].reshape(2*grid)                       # fill bulk with energy of nearest interface

    if options.ndimage:
      periodic_diffusedEnergy = ndimage.gaussian_filter(
                                np.where(ndimage.morphology.binary_dilation(periodic_interfaceEnergy > 0.,
                                                                            structure = struc,
                                                                            iterations = int(round(options.d*2.))-1,   # fat boundary
                                                                           ),
                                         periodic_bulkEnergy,                                      # ...and zero everywhere else
                                         0.),
                                sigma = options.d)
    else:
      diffusedEnergy = np.fft.irfftn(np.fft.rfftn(
                       np.where(
                         ndimage.morphology.binary_dilation(interfaceEnergy > 0.,
                                                            structure = struc,
                                                            iterations = int(round(options.d*2.))-1),# fat boundary
                         periodic_bulkEnergy[grid[0]//2:-grid[0]//2,                                 # retain filled energy on fat boundary...
                                             grid[1]//2:-grid[1]//2,
                                             grid[2]//2:-grid[2]//2],                               # ...and zero everywhere else
                         0.)).astype(np.complex64) *
                         gauss).astype(np.float32)

      periodic_diffusedEnergy = np.tile(diffusedEnergy,(3,3,3))[grid[0]//2:-grid[0]//2,
                                                                grid[1]//2:-grid[1]//2,
                                                                grid[2]//2:-grid[2]//2]             # periodically extend the smoothed bulk energy


    # transform voxels close to interface region
    index = ndimage.morphology.distance_transform_edt(periodic_diffusedEnergy >= 0.95*np.amax(periodic_diffusedEnergy),
                                                      return_distances = False,
                                                      return_indices = True)                        # want index of closest bulk grain

    periodic_material = np.tile(material,(3,3,3))[grid[0]//2:-grid[0]//2,
                                                  grid[1]//2:-grid[1]//2,
                                                  grid[2]//2:-grid[2]//2]                           # periodically extend the geometry

    material = periodic_material[index[0],
                                index[1],
                                index[2]].reshape(2*grid)[grid[0]//2:-grid[0]//2,
                                                          grid[1]//2:-grid[1]//2,
                                                          grid[2]//2:-grid[2]//2]                   # extent grains into interface region

    # replace immutable materials with closest mutable ones
    index = ndimage.morphology.distance_transform_edt(np.in1d(material,options.immutable).reshape(grid),
                                                      return_distances = False,
                                                      return_indices = True)
    material = material[index[0],
                          index[1],
                          index[2]]

    immutable = np.zeros(material.shape, dtype=np.bool)
    # find locations where immutable materials have been in original structure
    for micro in options.immutable:
      immutable += material_original == micro

    # undo any changes involving immutable materials
    material = np.where(immutable, material_original,material)

  damask.Grid(material = material[0:grid_original[0],0:grid_original[1],0:grid_original[2]],
              size      = geom.size,
              origin    = geom.origin,
              comments  = geom.comments + [scriptID + ' ' + ' '.join(sys.argv[1:])],
             )\
        .save(sys.stdout if name is None else name)
