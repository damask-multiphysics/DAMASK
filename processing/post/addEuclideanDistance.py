#!/usr/bin/env python2.7
# -*- coding: UTF-8 no BOM -*-

import os,sys,itertools
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
  for side in range(3):
    for p[(side+2)%3] in range(padded.shape[(side+2)%3]):
      for p[(side+1)%3] in range(padded.shape[(side+1)%3]):
        for p[side%3] in range(rimdim[side%3]):
          spot = (p-rimdim)%size
          padded[p[0],p[1],p[2]] = array[spot[0],spot[1],spot[2]]
        for p[side%3] in range(rimdim[side%3]+size[side%3],size[side%3]+2*rimdim[side%3]):
          spot = (p-rimdim)%size
          padded[p[0],p[1],p[2]] = array[spot[0],spot[1],spot[2]]
  return padded

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

features = [
            {'aliens': 1, 'names': ['boundary','biplane'],},
            {'aliens': 2, 'names': ['tripleline',],},
            {'aliens': 3, 'names': ['quadruplepoint',],}
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
Add column(s) containing Euclidean distance to grain structural features: boundaries, triple lines, and quadruple points.

""", version = scriptID)

parser.add_option('-p',
                  '--pos', '--position',
                  dest = 'pos', metavar = 'string',
                  help = 'label of coordinates [%default]')
parser.add_option('-i',
                  '--id', '--identifier',
                  dest = 'id', metavar = 'string',
                  help='label of grain identifier [%default]')
parser.add_option('-t',
                  '--type',
                  dest = 'type', action = 'extend', metavar = '<string LIST>',
                  help = 'feature type {{{}}} '.format(', '.join(map(lambda x:'/'.join(x['names']),features))) )
parser.add_option('-n',
                  '--neighborhood',
                  dest = 'neighborhood', choices = neighborhoods.keys(), metavar = 'string',
                  help = 'neighborhood type [neumann] {{{}}}'.format(', '.join(neighborhoods.keys())))
parser.add_option('-s',
                  '--scale',
                  dest = 'scale', type = 'float', metavar = 'float',
                  help = 'voxel size [%default]')

parser.set_defaults(pos = 'pos',
                    id = 'texture',
                    neighborhood = 'neumann',
                    scale = 1.0,
                   )

(options,filenames) = parser.parse_args()

if options.type is None:
  parser.error('no feature type selected.')
if not set(options.type).issubset(set(list(itertools.chain(*map(lambda x: x['names'],features))))):
  parser.error('type must be chosen from (%s).'%(', '.join(map(lambda x:'|'.join(x['names']),features))) )
if 'biplane' in options.type and 'boundary' in options.type:
  parser.error('only one from aliases "biplane" and "boundary" possible.')

feature_list = []
for i,feature in enumerate(features):
  for name in feature['names']:
    for myType in options.type:
      if name.startswith(myType):
        feature_list.append(i)                                                                      # remember valid features
        break

# --- loop over input files -------------------------------------------------------------------------

if filenames == []: filenames = [None]

for name in filenames:
  try:    table = damask.ASCIItable(name = name, buffered = False)
  except: continue
  damask.util.report(scriptName,name)

# ------------------------------------------ read header ------------------------------------------

  table.head_read()

# ------------------------------------------ sanity checks ----------------------------------------

  errors  = []
  remarks = []
  column = {}
  
  coordDim = table.label_dimension(options.pos)
  if not 3 >= coordDim >= 1:
    errors.append('coordinates "{}" need to have one, two, or three dimensions.'.format(options.pos))
  else: coordCol = table.label_index(options.pos)
  
  if table.label_dimension(options.id) != 1: errors.append('grain identifier {} not found.'.format(options.id))
  else: idCol = table.label_index(options.id)

  if remarks != []:
    damask.util.croak(remarks)
    remarks = []
  if errors  != []:
    damask.util.croak(errors)
    table.close(dismiss = True)
    continue

# ------------------------------------------ assemble header ---------------------------------------

  table.info_append(scriptID + '\t' + ' '.join(sys.argv[1:]))
  for feature in feature_list:
    table.labels_append('ED_{}({})'.format(features[feature]['names'][0],options.id))               # extend ASCII header with new labels
  table.head_write()

# --------------- figure out size and grid ---------------------------------------------------------

  table.data_readArray()

  coords = [np.unique(table.data[:,coordCol+i]) for i in range(coordDim)]
  mincorner = np.array(map(min,coords))
  maxcorner = np.array(map(max,coords))
  grid   = np.array(map(len,coords)+[1]*(3-len(coords)),'i')

  N = grid.prod()

  if N != len(table.data): errors.append('data count {} does not match grid {}.'.format(N,'x'.join(map(str,grid))))
  else:                   remarks.append('grid: {}x{}x{}'.format(*grid))
  if remarks != []: damask.util.croak(remarks)
  if errors  != []:
    damask.util.croak(errors)
    table.close(dismiss = True)
    continue

# ------------------------------------------ process value field -----------------------------------

  stack = [table.data]

  neighborhood = neighborhoods[options.neighborhood]
  diffToNeighbor = np.empty(list(grid+2)+[len(neighborhood)],'i')
  microstructure = periodic_3Dpad(table.data[:,idCol].astype('i').reshape(grid,order='F'))

  for i,p in enumerate(neighborhood):
    stencil = np.zeros((3,3,3),'i')
    stencil[1,1,1] = -1
    stencil[p[0]+1,
            p[1]+1,
            p[2]+1] = 1
    diffToNeighbor[:,:,:,i] = ndimage.convolve(microstructure,stencil)                              # compare ID at each point...
                                                                                                    # ...to every one in the specified neighborhood
                                                                                                    # for same IDs at both locations ==> 0
  
  diffToNeighbor = np.sort(diffToNeighbor)                                                          # sort diff such that number of changes in diff (steps)...
                                                                                                    # ...reflects number of unique neighbors
  uniques = np.where(diffToNeighbor[1:-1,1:-1,1:-1,0] != 0, 1,0)                                    # initialize unique value counter (exclude myself [= 0])

  for i in range(1,len(neighborhood)):                                                              # check remaining points in neighborhood
    uniques += np.where(np.logical_and(
                         diffToNeighbor[1:-1,1:-1,1:-1,i] != 0,                                     # not myself?
                         diffToNeighbor[1:-1,1:-1,1:-1,i] != diffToNeighbor[1:-1,1:-1,1:-1,i-1],
                         ),                                                                         # flip of ID difference detected?
                        1,0)                                                                        # count that flip

  distance = np.ones((len(feature_list),grid[0],grid[1],grid[2]),'d')

  for i,feature_id in enumerate(feature_list):
    distance[i,:,:,:] = np.where(uniques >= features[feature_id]['aliens'],0.0,1.0)                 # seed with 0.0 when enough unique neighbor IDs are present
    distance[i,:,:,:] = ndimage.morphology.distance_transform_edt(distance[i,:,:,:])*[options.scale]*3

  distance = distance.reshape([len(feature_list),grid.prod(),1],order='F')
  for i in range(len(feature_list)):
    stack.append(distance[i,:])

# ------------------------------------------ output result -----------------------------------------

  if len(stack) > 1: table.data = np.hstack(tuple(stack))
  table.data_writeArray('%.12g')

# ------------------------------------------ output finalization -----------------------------------

  table.close()                                                                                     # close input ASCII table (works for stdin)
