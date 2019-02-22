#!/usr/bin/env python3
# -*- coding: UTF-8 no BOM -*-

import os,sys,math
import numpy as np
from optparse import OptionParser
from scipy import ndimage
import damask

scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])

#--------------------------------------------------------------------------------------------------
#                                MAIN
#--------------------------------------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog [option(s)] [geomfile(s)]', description = """
Smoothens out interface roughness by simulated curvature flow.
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
                  help = 'list of immutable microstructure indices')
parser.add_option('-r', '--renumber',
                  dest = 'renumber', action='store_true',
                  help = 'output consecutive microstructure indices')
parser.add_option('--ndimage',
                  dest = 'ndimage', action='store_true',
                  help = 'use ndimage.gaussian_filter in lieu of explicit FFT')

parser.set_defaults(d = 1,
                    N = 1,
                   immutable = [],
                   renumber = False,
                   ndimage = False,
                   )

(options, filenames) = parser.parse_args()

options.immutable = list(map(int,options.immutable))

getInterfaceEnergy = lambda A,B: np.float32((A*B != 0)*(A != B)*1.0)                               # 1.0 if A & B are distinct & nonzero, 0.0 otherwise
struc = ndimage.generate_binary_structure(3,1)                                                     # 3D von Neumann neighborhood

# --- loop over input files -----------------------------------------------------------------------

if filenames == []: filenames = [None]

for name in filenames:
  try:    table = damask.ASCIItable(name = name,
                                    buffered = False,
                                    labeled = False)
  except: continue
  damask.util.report(scriptName,name)

# --- interpret header ----------------------------------------------------------------------------

  table.head_read()
  info,extra_header = table.head_getGeom()
  
  damask.util.croak(['grid     a b c:  {}'.format(' x '.join(list(map(str,info['grid'])))),
                     'size     x y z:  {}'.format(' x '.join(list(map(str,info['size'])))),
                     'origin   x y z:  {}'.format(' : '.join(list(map(str,info['origin'])))),
                     'homogenization:  {}'.format(info['homogenization']),
                     'microstructures: {}'.format(info['microstructures']),
                    ])

  errors = []
  if np.any(info['grid'] < 1):    errors.append('invalid grid a b c.')
  if np.any(info['size'] <= 0.0): errors.append('invalid size x y z.')
  if errors != []:
    damask.util.croak(errors)
    table.close(dismiss = True)
    continue

# --- read data -----------------------------------------------------------------------------------
  microstructure = np.tile(table.microstructure_read(info['grid']).reshape(info['grid'],order='F'),
                           np.where(info['grid'] == 1, 2,1))                                       # make one copy along dimensions with grid == 1
  grid = np.array(microstructure.shape)

# --- initialize support data ---------------------------------------------------------------------

# store a copy the initial microstructure to find locations of immutable indices
  microstructure_original = np.copy(microstructure)                                                 

  if not options.ndimage:
    X,Y,Z = np.mgrid[0:grid[0],0:grid[1],0:grid[2]]
  
    # Calculates gaussian weights for simulating 3d diffusion
    gauss = np.exp(-(X*X + Y*Y + Z*Z)/(2.0*options.d*options.d),dtype=np.float32) \
            /np.power(2.0*np.pi*options.d*options.d,(3.0 - np.count_nonzero(info['grid'] == 1))/2.,dtype=np.float32)
          
    gauss[:,:,:grid[2]//2:-1] = gauss[:,:,1:(grid[2]+1)//2]     # trying to cope with uneven (odd) grid size
    gauss[:,:grid[1]//2:-1,:] = gauss[:,1:(grid[1]+1)//2,:]
    gauss[:grid[0]//2:-1,:,:] = gauss[1:(grid[0]+1)//2,:,:]
    gauss = np.fft.rfftn(gauss).astype(np.complex64)

  for smoothIter in range(options.N):

    interfaceEnergy = np.zeros(microstructure.shape,dtype=np.float32)
    for i in (-1,0,1):
      for j in (-1,0,1):
        for k in (-1,0,1):
          # assign interfacial energy to all voxels that have a differing neighbor (in Moore neighborhood)
          interfaceEnergy = np.maximum(interfaceEnergy,
                                       getInterfaceEnergy(microstructure,np.roll(np.roll(np.roll(
                                                          microstructure,i,axis=0), j,axis=1), k,axis=2)))

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
                                             grid[2]//2:-grid[2]//2],                                # ...and zero everywhere else
                         0.)).astype(np.complex64) *
                         gauss).astype(np.float32)

      periodic_diffusedEnergy = np.tile(diffusedEnergy,(3,3,3))[grid[0]//2:-grid[0]//2,
                                                                grid[1]//2:-grid[1]//2,
                                                                grid[2]//2:-grid[2]//2]              # periodically extend the smoothed bulk energy


    # transform voxels close to interface region
    index = ndimage.morphology.distance_transform_edt(periodic_diffusedEnergy >= 0.95*np.amax(periodic_diffusedEnergy),
                                                      return_distances = False,
                                                      return_indices = True)                       # want index of closest bulk grain

    periodic_microstructure = np.tile(microstructure,(3,3,3))[grid[0]//2:-grid[0]//2,
                                                              grid[1]//2:-grid[1]//2,
                                                              grid[2]//2:-grid[2]//2]                # periodically extend the microstructure

    microstructure = periodic_microstructure[index[0],
                                             index[1],
                                             index[2]].reshape(2*grid)[grid[0]//2:-grid[0]//2,
                                                                       grid[1]//2:-grid[1]//2,
                                                                       grid[2]//2:-grid[2]//2]       # extent grains into interface region

    # replace immutable microstructures with closest mutable ones
    index = ndimage.morphology.distance_transform_edt(np.in1d(microstructure,options.immutable).reshape(grid),
                                                      return_distances = False,
                                                      return_indices = True)
    microstructure = microstructure[index[0],
                                    index[1],
                                    index[2]]

    immutable = np.zeros(microstructure.shape, dtype=np.bool)
    # find locations where immutable microstructures have been in original structure
    for micro in options.immutable:
      immutable += microstructure_original == micro

    # undo any changes involving immutable microstructures
    microstructure = np.where(immutable, microstructure_original,microstructure)

# --- renumber to sequence 1...Ngrains if requested -----------------------------------------------
#  http://stackoverflow.com/questions/10741346/np-frequency-counts-for-unique-values-in-an-array

  if options.renumber:
    newID = 0
    for microstructureID,count in enumerate(np.bincount(microstructure.flatten())):
      if count != 0:
        newID += 1
        microstructure = np.where(microstructure == microstructureID, newID, microstructure)

  newInfo = {'microstructures': 0,}
  newInfo['microstructures'] = microstructure.max()

# --- report --------------------------------------------------------------------------------------

  remarks = []
  if newInfo['microstructures'] != info['microstructures']:
    remarks.append('--> microstructures: {}'.format(newInfo['microstructures']))
  if remarks != []: damask.util.croak(remarks)

# --- write header --------------------------------------------------------------------------------

  table.labels_clear()
  table.info_clear()
  table.info_append(extra_header+[
    scriptID + ' ' + ' '.join(sys.argv[1:]),
    "grid\ta {grid[0]}\tb {grid[1]}\tc {grid[2]}".format(grid=info['grid']),
    "size\tx {size[0]}\ty {size[1]}\tz {size[2]}".format(size=info['size']),
    "origin\tx {origin[0]}\ty {origin[1]}\tz {origin[2]}".format(origin=info['origin']),
    "homogenization\t{homog}".format(homog=info['homogenization']),
    "microstructures\t{microstructures}".format(microstructures=newInfo['microstructures']),
    ])
  table.head_write()

# --- write microstructure information ------------------------------------------------------------

  formatwidth = int(math.floor(math.log10(microstructure.max())+1))
  table.data = microstructure[::1 if info['grid'][0]>1 else 2,
                              ::1 if info['grid'][1]>1 else 2,
                              ::1 if info['grid'][2]>1 else 2,].\
                              reshape((info['grid'][0],info['grid'][1]*info['grid'][2]),order='F').transpose()
  table.data_writeArray('%{}i'.format(formatwidth),delimiter = ' ')

# --- output finalization --------------------------------------------------------------------------

  table.close()     
  