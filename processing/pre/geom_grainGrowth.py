#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,string,re,math,itertools
import numpy as np
from optparse import OptionParser
from scipy import ndimage
from multiprocessing import Pool
import damask

scriptID = '$Id$'
scriptName = scriptID.split()[1]

#--------------------------------------------------------------------------------------------------
#                                MAIN
#--------------------------------------------------------------------------------------------------
synonyms = {
        'grid':   ['resolution'],
        'size':   ['dimension'],
          }
identifiers = {
        'grid':   ['a','b','c'],
        'size':   ['x','y','z'],
        'origin': ['x','y','z'],
          }
mappings = {
        'grid':            lambda x: int(x),
        'size':            lambda x: float(x),
        'origin':          lambda x: float(x),
        'homogenization':  lambda x: int(x),
        'microstructures': lambda x: int(x),
          }

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Smoothens out interface roughness by simulated curvature flow.
This is achieved by the diffusion of each initially sharply bounded grain volume within the periodic domain
up to a given distance 'd' voxels.
The final geometry is assembled by selecting at each voxel that grain index for which the concentration remains largest.
""" + string.replace(scriptID,'\n','\\n')
)

parser.add_option('-d', '--distance', dest='d', type='int', metavar='int',
                  help='diffusion distance in voxels [%default]')
parser.add_option('-N', '--smooth', dest='N', type='int', metavar='int',
                 help='N for curvature flow [%default]')
parser.add_option('-r', '--renumber', dest='renumber', action='store_true',
                  help='renumber microstructure indices from 1...N [%default]')
parser.add_option('-i', '--immutable', action='extend', dest='immutable', type='string', metavar = '<LIST>',
                  help='list of immutable microstructures')

parser.set_defaults(d = 1)
parser.set_defaults(N = 1)
parser.set_defaults(renumber = False)
parser.set_defaults(immutable = [])

(options, filenames) = parser.parse_args()

options.immutable = map(int,options.immutable)

#--- setup file handles --------------------------------------------------------------------------   
files = []
if filenames == []:
 files.append({'name':'STDIN',
               'input':sys.stdin,
               'output':sys.stdout,
               'croak':sys.stderr,
              })
else:
 for name in filenames:
   if os.path.exists(name):
     files.append({'name':name,
                   'input':open(name),
                   'output':open(name+'_tmp','w'),
                   'croak':sys.stdout,
                   })

#--- loop over input files ------------------------------------------------------------------------
for file in files:
  if file['name'] != 'STDIN': file['croak'].write('\033[1m'+scriptName+'\033[0m: '+file['name']+'\n')
  else: file['croak'].write('\033[1m'+scriptName+'\033[0m\n')

  theTable = damask.ASCIItable(file['input'],file['output'],labels = False,buffered = False)
  theTable.head_read()

#--- interpret header ----------------------------------------------------------------------------
  info = {
          'grid':   np.zeros(3,'i'),
          'size':   np.zeros(3,'d'),
          'origin': np.zeros(3,'d'),
          'homogenization':  0,
          'microstructures': 0,
         }
  newInfo = {
          'microstructures': 0,
         }
  extra_header = []

  for header in theTable.info:
    headitems = map(str.lower,header.split())
    if len(headitems) == 0: continue
    for synonym,alternatives in synonyms.iteritems():
      if headitems[0] in alternatives: headitems[0] = synonym
    if headitems[0] in mappings.keys():
      if headitems[0] in identifiers.keys():
        for i in xrange(len(identifiers[headitems[0]])):
          info[headitems[0]][i] = \
            mappings[headitems[0]](headitems[headitems.index(identifiers[headitems[0]][i])+1])
      else:
        info[headitems[0]] = mappings[headitems[0]](headitems[1])
    else:
      extra_header.append(header)

  file['croak'].write('grid     a b c:  %s\n'%(' x '.join(map(str,info['grid']))) + \
                      'size     x y z:  %s\n'%(' x '.join(map(str,info['size']))) + \
                      'origin   x y z:  %s\n'%(' : '.join(map(str,info['origin']))) + \
                      'homogenization:  %i\n'%info['homogenization'] + \
                      'microstructures: %i\n'%info['microstructures'])

  if np.any(info['grid'] < 1):
    file['croak'].write('invalid grid a b c.\n')
    continue
  if np.any(info['size'] <= 0.0):
    file['croak'].write('invalid size x y z.\n')
    continue

#--- read data ------------------------------------------------------------------------------------
  microstructure = np.zeros(np.prod(info['grid']),'i')                                                # 2D structures do not work
  i = 0

  while theTable.data_read():                                  # read next data line of ASCII table
    items = theTable.data
    if len(items) > 2:
      if   items[1].lower() == 'of': items = [int(items[2])]*int(items[0])
      elif items[1].lower() == 'to': items = xrange(int(items[0]),1+int(items[2]))
      else:                          items = map(int,items)
    else:                            items = map(int,items)

    s = len(items)
    microstructure[i:i+s] = items
    i += s

#--- reshape, if 2D make copy ---------------------------------------------------------------------

  microstructure = np.tile(microstructure.reshape(info['grid'],order='F'),
                              np.where(info['grid'] == 1, 2,1))                                                                     # make one copy along dimensions with grid == 1
  grid = np.array(microstructure.shape)
  
#--- initialize support data -----------------------------------------------------------------------  

  periodic_microstructure = np.tile(microstructure,(3,3,3))[grid[0]/2:-grid[0]/2,
                                                               grid[1]/2:-grid[1]/2,
                                                               grid[2]/2:-grid[2]/2]                                                # periodically extend the microstructure
  microstructure_original = np.copy(microstructure)                                                                                 # store a copy the initial microstructure to find locations of immutable indices

  X,Y,Z = np.mgrid[0:grid[0],0:grid[1],0:grid[2]]
  gauss = np.exp(-(X*X + Y*Y + Z*Z)/(2.0*options.d*options.d))/math.pow(2.0*np.pi*options.d*options.d,1.5)
  gauss[:,:,grid[2]/2::] = gauss[:,:,round(grid[2]/2.)-1::-1]     # trying to cope with uneven (odd) grid size
  gauss[:,grid[1]/2::,:] = gauss[:,round(grid[1]/2.)-1::-1,:]
  gauss[grid[0]/2::,:,:] = gauss[round(grid[0]/2.)-1::-1,:,:]
  gauss = np.fft.rfftn(gauss)
  
  interfacialEnergy = lambda A,B: (A*B != 0)*(A != B)*1.0
  struc = ndimage.generate_binary_structure(3,1)                                                                                    # 3D von Neumann neighborhood


  for smoothIter in xrange(options.N):
    boundary = np.zeros(microstructure.shape)
    for i in (-1,0,1):
      for j in (-1,0,1):
        for k in (-1,0,1):
          interfaceEnergy = np.maximum(boundary,
                                          interfacialEnergy(microstructure,np.roll(np.roll(np.roll(
                                                            microstructure,i,axis=0), j,axis=1), k,axis=2)))                        # assign interfacial energy to all voxels that have a differing neighbor (in Moore neighborhood)
    periodic_interfaceEnergy = np.tile(interfaceEnergy,(3,3,3))[grid[0]/2:-grid[0]/2,
                                                                   grid[1]/2:-grid[1]/2,
                                                                   grid[2]/2:-grid[2]/2]                                            # periodically extend interfacial energy array by half a grid size in positive and negative directions
    index = ndimage.morphology.distance_transform_edt(periodic_interfaceEnergy == 0.,                                               # transform bulk volume (i.e. where interfacial energy is zero)
                                                      return_distances = False,
                                                      return_indices = True)                                                        # want array index of nearest voxel on periodically extended boundary
#    boundaryExt = boundaryExt[index[0].flatten(),index[1].flatten(),index[2].flatten()].reshape(boundaryExt.shape)                 # fill bulk with energy of nearest interface | question PE: what "flatten" for?
    periodic_bulkEnergy = periodic_interfaceEnergy[index[0],
                                                   index[1],
                                                   index[2]].reshape(2*grid)                                                        # fill bulk with energy of nearest interface
    diffusedEnergy = np.fft.irfftn(np.fft.rfftn(np.where(ndimage.morphology.binary_dilation(interfaceEnergy > 0.,
                                                                                                     structure = struc,
                                                                                                     iterations = options.d/2 + 1), # fat boundary | question PE: why 2d - 1? I would argue for d/2 + 1 !!
                                                            periodic_bulkEnergy[grid[0]/2:-grid[0]/2,                               # retain filled energy on fat boundary...
                                                                                grid[1]/2:-grid[1]/2,
                                                                                grid[2]/2:-grid[2]/2],                              # ...and zero everywhere else
                                                            0.)\
                                                )*gauss)
    periodic_diffusedEnergy = np.tile(diffusedEnergy,(3,3,3))[grid[0]/2:-grid[0]/2,
                                                                 grid[1]/2:-grid[1]/2,
                                                                 grid[2]/2:-grid[2]/2]                                              # periodically extend the smoothed bulk energy
    index = ndimage.morphology.distance_transform_edt(periodic_diffusedEnergy >= 0.5,                                               # transform voxels close to interface region | question PE: what motivates 1/2 (could be any small number, or)?
                                                      return_distances = False,
                                                      return_indices = True)                                                        # want index of closest bulk grain
    microstructure = periodic_microstructure[index[0],
                                             index[1],
                                             index[2]].reshape(2*grid)[grid[0]/2:-grid[0]/2,
                                                                       grid[1]/2:-grid[1]/2,
                                                                       grid[2]/2:-grid[2]/2]                                        # extent grains into interface region

    immutable = np.zeros(microstructure.shape, dtype=bool)
    for micro in options.immutable:
      immutable += np.logical_or(microstructure == micro, microstructure_original == micro)                                         # find locations where immutable microstructures have been or are now

    microstructure = np.where(immutable, microstructure_original,microstructure)                                                    # undo any changes involving immutable microstructures

# --- renumber to sequence 1...Ngrains if requested ------------------------------------------------
#  http://stackoverflow.com/questions/10741346/np-frequency-counts-for-unique-values-in-an-array  

  if options.renumber:
    newID = 0
    for microstructureID,count in enumerate(np.bincount(microstructure.flatten())):
      if count != 0:
        newID += 1
        microstructure = np.where(microstructure == microstructureID, newID, microstructure)

# --- assemble header -----------------------------------------------------------------------------
  newInfo['microstructures'] = microstructure[0:info['grid'][0],0:info['grid'][1],0:info['grid'][2]].max()

#--- report ---------------------------------------------------------------------------------------
  if (newInfo['microstructures'] != info['microstructures']):
    file['croak'].write('--> microstructures: %i\n'%newInfo['microstructures'])

#--- write header ---------------------------------------------------------------------------------
  theTable.labels_clear()
  theTable.info_clear()
  theTable.info_append(extra_header+[
    scriptID+ ' ' + ' '.join(sys.argv[1:]),
    "grid\ta %i\tb %i\tc %i"%(info['grid'][0],info['grid'][1],info['grid'][2],),
    "size\tx %f\ty %f\tz %f"%(info['size'][0],info['size'][1],info['size'][2],),
    "origin\tx %f\ty %f\tz %f"%(info['origin'][0],info['origin'][1],info['origin'][2],),
    "homogenization\t%i"%info['homogenization'],
    "microstructures\t%i"%(newInfo['microstructures']),
    ])
  theTable.head_write()
  
# --- write microstructure information ------------------------------------------------------------
  formatwidth = int(math.floor(math.log10(microstructure.max())+1))
  theTable.data = microstructure[0:info['grid'][0],0:info['grid'][1],0:info['grid'][2]].reshape(np.prod(info['grid']),order='F').transpose() # question PE: this assumes that only the Z dimension can be 1!
  theTable.data_writeArray('%%%ii'%(formatwidth),delimiter=' ')
    
#--- output finalization --------------------------------------------------------------------------
  if file['name'] != 'STDIN':
    theTable.__IO__['in'].close()
    theTable.__IO__['out'].close()
    os.rename(file['name']+'_tmp',file['name'])
