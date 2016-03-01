#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,math,itertools
import numpy as np
from scipy import ndimage
from optparse import OptionParser
import damask

scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])

def periodic_3Dpad(array, rimdim=(1,1,1)):

  rimdim = np.array(rimdim,'i')
  size = np.array(array.shape,'i')
  padded = np.empty(size+2*rimdim,array.dtype)
  padded[rimdim[0]:rimdim[0]+size[0],
         rimdim[1]:rimdim[1]+size[1],
         rimdim[2]:rimdim[2]+size[2]] = array

  p = np.zeros(3,'i')
  for side in xrange(3):
    for p[(side+2)%3] in xrange(padded.shape[(side+2)%3]):
      for p[(side+1)%3] in xrange(padded.shape[(side+1)%3]):
        for p[side%3] in xrange(rimdim[side%3]):
          spot = (p-rimdim)%size
          padded[p[0],p[1],p[2]] = array[spot[0],spot[1],spot[2]]
        for p[side%3] in xrange(rimdim[side%3]+size[side%3],size[side%3]+2*rimdim[side%3]):
          spot = (p-rimdim)%size
          padded[p[0],p[1],p[2]] = array[spot[0],spot[1],spot[2]]
  return padded

#--------------------------------------------------------------------------------------------------
#                                MAIN
#--------------------------------------------------------------------------------------------------

features = [
            {'aliens': 1, 'alias': ['boundary','biplane'],},
            {'aliens': 2, 'alias': ['tripleline',],},
            {'aliens': 3, 'alias': ['quadruplepoint',],}
           ]

neighborhoods = {
                  'neumann':np.array([
                                        [-1, 0, 0],
                                        [ 1, 0, 0],
                                        [ 0,-1, 0],
                                        [ 0, 1, 0],
                                        [ 0, 0,-1],
                                        [ 0, 0, 1],
                                      ]),
                  'moore':np.array([
                                        [-1,-1,-1],
                                        [ 0,-1,-1],
                                        [ 1,-1,-1],
                                        [-1, 0,-1],
                                        [ 0, 0,-1],
                                        [ 1, 0,-1],
                                        [-1, 1,-1],
                                        [ 0, 1,-1],
                                        [ 1, 1,-1],
#
                                        [-1,-1, 0],
                                        [ 0,-1, 0],
                                        [ 1,-1, 0],
                                        [-1, 0, 0],
#
                                        [ 1, 0, 0],
                                        [-1, 1, 0],
                                        [ 0, 1, 0],
                                        [ 1, 1, 0],
#
                                        [-1,-1, 1],
                                        [ 0,-1, 1],
                                        [ 1,-1, 1],
                                        [-1, 0, 1],
                                        [ 0, 0, 1],
                                        [ 1, 0, 1],
                                        [-1, 1, 1],
                                        [ 0, 1, 1],
                                        [ 1, 1, 1],
                                      ])
                }

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Produce geom files containing Euclidean distance to grain structural features:
boundaries, triple lines, and quadruple points.

""", version = scriptID)

parser.add_option('-t','--type',
                  dest = 'type',
                  action = 'extend', metavar = '<string LIST>',
                  help = 'feature type {%s} '%(', '.join(map(lambda x:'|'.join(x['alias']),features))) )
parser.add_option('-n','--neighborhood',
                  dest = 'neighborhood',
                  choices = neighborhoods.keys(), metavar = 'string',
                  help = 'type of neighborhood {%s} [neumann]'%(', '.join(neighborhoods.keys())))
parser.add_option('-s', '--scale',
                  dest = 'scale',
                  type = 'float', metavar = 'float',
                  help = 'voxel size [%default]')

parser.set_defaults(type = [],
                    neighborhood = 'neumann',
                    scale = 1.0,
                   )

(options,filenames) = parser.parse_args()

if len(options.type) == 0 or \
   not set(options.type).issubset(set(list(itertools.chain(*map(lambda x: x['alias'],features))))):
  parser.error('sleect feature type from (%s).'%(', '.join(map(lambda x:'|'.join(x['alias']),features))) )
if 'biplane' in options.type and 'boundary' in options.type:
  parser.error("only one alias out 'biplane' and 'boundary' required")
  
feature_list = []
for i,feature in enumerate(features):
  for name in feature['alias']:
    for myType in options.type:
      if name.startswith(myType):
        feature_list.append(i)                                                                      # remember selected features
        break

# --- loop over input files -------------------------------------------------------------------------

if filenames == []: filenames = [None]

for name in filenames:
  try:
    table = damask.ASCIItable(name = name,
                              buffered = False, labeled = False, readonly = True)
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

  microstructure = table.microstructure_read(info['grid']).reshape(info['grid'],order='F')          # read microstructure

  table.close()
  
  neighborhood = neighborhoods[options.neighborhood]
  convoluted = np.empty([len(neighborhood)]+list(info['grid']+2),'i')
  structure = periodic_3Dpad(microstructure)
  
  for i,p in enumerate(neighborhood):
    stencil = np.zeros((3,3,3),'i')
    stencil[1,1,1] = -1
    stencil[p[0]+1,
            p[1]+1,
            p[2]+1] = 1
    convoluted[i,:,:,:] = ndimage.convolve(structure,stencil)
  
  convoluted = np.sort(convoluted,axis = 0)
  uniques = np.where(convoluted[0,1:-1,1:-1,1:-1] != 0, 1,0)                                        # initialize unique value counter (exclude myself [= 0])

  for i in xrange(1,len(neighborhood)):                                                             # check remaining points in neighborhood
    uniques += np.where(np.logical_and(
                           convoluted[i,1:-1,1:-1,1:-1] != convoluted[i-1,1:-1,1:-1,1:-1],          # flip of ID difference detected?
                           convoluted[i,1:-1,1:-1,1:-1] != 0),                                      # not myself?
                           1,0)                                                                     # count flip

  for feature in feature_list:
    try:
      table = damask.ASCIItable(outname = features[feature]['alias'][0]+'_'+name if name else name,
                                buffered = False, labeled = False)
    except: continue

    damask.util.croak(features[feature]['alias'][0])
      
    distance = np.where(uniques >= features[feature]['aliens'],0.0,1.0)                             # seed with 0.0 when enough unique neighbor IDs are present
    distance = ndimage.morphology.distance_transform_edt(distance)*[options.scale]*3

    info['microstructures'] = int(math.ceil(distance.max()))

#--- write header ---------------------------------------------------------------------------------

    table.info_clear()
    table.info_append(extra_header+[
      scriptID + ' ' + ' '.join(sys.argv[1:]),
      "grid\ta {grid[0]}\tb {grid[1]}\tc {grid[2]}".format(grid=info['grid']),
      "size\tx {size[0]}\ty {size[1]}\tz {size[2]}".format(size=info['size']),
      "origin\tx {origin[0]}\ty {origin[1]}\tz {origin[2]}".format(origin=info['origin']),
      "homogenization\t{homog}".format(homog=info['homogenization']),
      "microstructures\t{microstructures}".format(microstructures=info['microstructures']),
      ])
    table.labels_clear()
    table.head_write()
    
# --- write microstructure information ------------------------------------------------------------

    formatwidth = int(math.floor(math.log10(distance.max())+1))
    table.data = distance.reshape((info['grid'][0],info['grid'][1]*info['grid'][2]),order='F').transpose()
    table.data_writeArray('%%%ii'%(formatwidth),delimiter=' ')
    
#--- output finalization --------------------------------------------------------------------------

    table.close()
