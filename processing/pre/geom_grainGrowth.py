#!/usr/bin/env python
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

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Smoothens out interface roughness by simulated curvature flow.
This is achieved by the diffusion of each initially sharply bounded grain volume within the periodic domain
up to a given distance 'd' voxels.
The final geometry is assembled by selecting at each voxel that grain index for which the concentration remains largest.

""", version = scriptID)

parser.add_option('-d', '--distance', dest='d', type='int', metavar='int',
                  help='diffusion distance in voxels [%default]')
parser.add_option('-N', '--smooth', dest='N', type='int', metavar='int',
                 help='N for curvature flow [%default]')
parser.add_option('-r', '--renumber', dest='renumber', action='store_true',
                  help='renumber microstructure indices from 1...N [%default]')
parser.add_option('-i', '--immutable', action='extend', dest='immutable', metavar = '<LIST>',
                  help='list of immutable microstructures')

parser.set_defaults(d = 1)
parser.set_defaults(N = 1)
parser.set_defaults(renumber = False)
parser.set_defaults(immutable = [])

(options, filenames) = parser.parse_args()

options.immutable = map(int,options.immutable)

# --- loop over input files -------------------------------------------------------------------------

if filenames == []: filenames = [None]

for name in filenames:
  try:
    table = damask.ASCIItable(name = name,
                              buffered = False, labeled = False)
  except: continue
  damask.util.report(scriptName,name)

# --- interpret header ----------------------------------------------------------------------------

  table.head_read()
  info,extra_header = table.head_getGeom()
  
  damask.util.croak(['grid     a b c:  %s'%(' x '.join(map(str,info['grid']))),
               'size     x y z:  %s'%(' x '.join(map(str,info['size']))),
               'origin   x y z:  %s'%(' : '.join(map(str,info['origin']))),
               'homogenization:  %i'%info['homogenization'],
               'microstructures: %i'%info['microstructures'],
              ])

  errors = []
  if np.any(info['grid'] < 1):    errors.append('invalid grid a b c.')
  if np.any(info['size'] <= 0.0): errors.append('invalid size x y z.')
  if errors != []:
    damask.util.croak(errors)
    table.close(dismiss = True)
    continue

# --- read data ------------------------------------------------------------------------------------
  microstructure = np.tile(np.array(table.microstructure_read(info['grid']),'i').reshape(info['grid'],order='F'),
                              np.where(info['grid'] == 1, 2,1))                                     # make one copy along dimensions with grid == 1
  grid = np.array(microstructure.shape)

#--- initialize support data -----------------------------------------------------------------------

  periodic_microstructure = np.tile(microstructure,(3,3,3))[grid[0]/2:-grid[0]/2,
                                                            grid[1]/2:-grid[1]/2,
                                                            grid[2]/2:-grid[2]/2]                   # periodically extend the microstructure
# store a copy the initial microstructure to find locations of immutable indices
  microstructure_original = np.copy(microstructure)                                                 

  X,Y,Z = np.mgrid[0:grid[0],0:grid[1],0:grid[2]]
  gauss = np.exp(-(X*X + Y*Y + Z*Z)/(2.0*options.d*options.d))/math.pow(2.0*np.pi*options.d*options.d,1.5)
  gauss[:,:,grid[2]/2::] = gauss[:,:,round(grid[2]/2.)-1::-1]     # trying to cope with uneven (odd) grid size
  gauss[:,grid[1]/2::,:] = gauss[:,round(grid[1]/2.)-1::-1,:]
  gauss[grid[0]/2::,:,:] = gauss[round(grid[0]/2.)-1::-1,:,:]
  gauss = np.fft.rfftn(gauss)

  interfacialEnergy = lambda A,B: (A*B != 0)*(A != B)*1.0
  struc = ndimage.generate_binary_structure(3,1)                                                                                 # 3D von Neumann neighborhood


  for smoothIter in xrange(options.N):
    boundary = np.zeros(microstructure.shape)
    for i in (-1,0,1):
      for j in (-1,0,1):
        for k in (-1,0,1):
            # assign interfacial energy to all voxels that have a differing neighbor (in Moore neighborhood)
          interfaceEnergy = np.maximum(boundary,
                                          interfacialEnergy(microstructure,np.roll(np.roll(np.roll(
                                                            microstructure,i,axis=0), j,axis=1), k,axis=2)))
    # periodically extend interfacial energy array by half a grid size in positive and negative directions
    periodic_interfaceEnergy = np.tile(interfaceEnergy,(3,3,3))[grid[0]/2:-grid[0]/2,
                                                                   grid[1]/2:-grid[1]/2,
                                                                   grid[2]/2:-grid[2]/2]
    # transform bulk volume (i.e. where interfacial energy is zero)
    index = ndimage.morphology.distance_transform_edt(periodic_interfaceEnergy == 0.,                                            
                                                      return_distances = False,
                                                      return_indices = True)
    # want array index of nearest voxel on periodically extended boundary
    periodic_bulkEnergy = periodic_interfaceEnergy[index[0],
                                                   index[1],
                                                   index[2]].reshape(2*grid)                        # fill bulk with energy of nearest interface
    diffusedEnergy = np.fft.irfftn(np.fft.rfftn(
                     np.where(
                       ndimage.morphology.binary_dilation(interfaceEnergy > 0.,
                                                          structure = struc,
                                                          iterations = options.d/2 + 1),            # fat boundary | PE: why 2d-1? I would argue for d/2 + 1
                       periodic_bulkEnergy[grid[0]/2:-grid[0]/2,                                    # retain filled energy on fat boundary...
                                           grid[1]/2:-grid[1]/2,
                                           grid[2]/2:-grid[2]/2],                                   # ...and zero everywhere else
                       0.))*gauss)
    periodic_diffusedEnergy = np.tile(diffusedEnergy,(3,3,3))[grid[0]/2:-grid[0]/2,
                                                                 grid[1]/2:-grid[1]/2,
                                                                 grid[2]/2:-grid[2]/2]              # periodically extend the smoothed bulk energy
    # transform voxels close to interface region | question PE: what motivates 1/2 (could be any small number, or)?
    index = ndimage.morphology.distance_transform_edt(periodic_diffusedEnergy >= 0.5,                                            
                                                      return_distances = False,
                                                      return_indices = True)                        # want index of closest bulk grain
    microstructure = periodic_microstructure[index[0],
                                             index[1],
                                             index[2]].reshape(2*grid)[grid[0]/2:-grid[0]/2,
                                                                       grid[1]/2:-grid[1]/2,
                                                                       grid[2]/2:-grid[2]/2]        # extent grains into interface region

    immutable = np.zeros(microstructure.shape, dtype=bool)
    # find locations where immutable microstructures have been or are now
    for micro in options.immutable:
      immutable += np.logical_or(microstructure == micro, microstructure_original == micro)         
    # undo any changes involving immutable microstructures
    microstructure = np.where(immutable, microstructure_original,microstructure)                    

# --- renumber to sequence 1...Ngrains if requested ------------------------------------------------
#  http://stackoverflow.com/questions/10741346/np-frequency-counts-for-unique-values-in-an-array

  if options.renumber:
    newID = 0
    for microstructureID,count in enumerate(np.bincount(microstructure.flatten())):
      if count != 0:
        newID += 1
        microstructure = np.where(microstructure == microstructureID, newID, microstructure)

  newInfo = {'microstructures': 0,}
  newInfo['microstructures'] = microstructure.max()

# --- report ---------------------------------------------------------------------------------------

  remarks = []
  if (newInfo['microstructures'] != info['microstructures']): remarks.append('--> microstructures: %i'%newInfo['microstructures'])
  if remarks != []: damask.util.croak(remarks)

# --- write header ---------------------------------------------------------------------------------

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
  table.data = microstructure.reshape((info['grid'][0],info['grid'][1]*info['grid'][2]),order='F').transpose()
  table.data_writeArray('%%%ii'%(formatwidth),delimiter = ' ')

# --- output finalization --------------------------------------------------------------------------

  table.close()     
