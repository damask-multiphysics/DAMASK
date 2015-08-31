#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,string,re,math,itertools
import numpy as np
from scipy import ndimage
from optparse import OptionParser
import damask

scriptID   = string.replace('$Id$','\n','\\n')
scriptName = os.path.splitext(scriptID.split()[1])[0]

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

parser.add_option('-c','--coordinates', dest='coords', metavar='string',
                  help='column heading for coordinates [%default]')
parser.add_option('-i','--identifier',  dest='id', metavar = 'string',
                  help='heading of column containing grain identifier [%default]')
parser.add_option('-t','--type',        dest = 'type', action = 'extend', metavar = '<string LIST>',
                  help = 'feature type {%s} '%(', '.join(map(lambda x:'/'.join(x['names']),features))) )
parser.add_option('-n','--neighborhood',dest='neighborhood', choices = neighborhoods.keys(), metavar = 'string',
                  help = 'type of neighborhood [neumann] {%s}'%(', '.join(neighborhoods.keys())))
parser.add_option('-s', '--scale',      dest = 'scale', type = 'float', metavar='float',
                  help = 'voxel size [%default]')
parser.set_defaults(coords = 'ipinitialcoord')
parser.set_defaults(id = 'texture')
parser.set_defaults(neighborhood = 'neumann')
parser.set_defaults(scale = 1.0)

(options,filenames) = parser.parse_args()

if options.type == None:
  parser.error('no feature type selected.')
if not set(options.type).issubset(set(list(itertools.chain(*map(lambda x: x['names'],features))))):
  parser.error('type must be chosen from (%s).'%(', '.join(map(lambda x:'|'.join(x['names']),features))) )
if 'biplane' in options.type and 'boundary' in options.type:
  parser.error("only one from aliases 'biplane' and 'boundary' possible.")

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
  try:
    table = damask.ASCIItable(name = name, buffered = False)
  except:
    continue
  table.croak(damask.util.emph(scriptName)+(': '+name if name else ''))

# ------------------------------------------ read header ------------------------------------------

  table.head_read()

# ------------------------------------------ sanity checks ----------------------------------------

  errors  = []
  remarks = []
  column = {}
  
  if table.label_dimension(options.coords) != 3: errors.append('coordinates {} are not a vector.'.format(options.coords))
  else: coordCol = table.label_index(options.coords)

  if table.label_dimension(options.id) != 1: errors.append('grain identifier {} not found.'.format(options.id))
  else: idCol = table.label_index(options.id)

  if remarks != []: table.croak(remarks)
  if errors  != []:
    table.croak(errors)
    table.close(dismiss = True)
    continue

# ------------------------------------------ assemble header ---------------------------------------

  table.info_append(scriptID + '\t' + ' '.join(sys.argv[1:]))
  for feature in feature_list:
    table.labels_append('ED_{}({})'.format(features[feature]['names'][0],options.id))               # extend ASCII header with new labels
  table.head_write()

# --------------- figure out size and grid ---------------------------------------------------------

  table.data_readArray()

  coords = [{},{},{}]
  for i in xrange(len(table.data)):  
    for j in xrange(3):
      coords[j][str(table.data[i,coordCol+j])] = True
  grid = np.array(map(len,coords),'i')
  size = grid/np.maximum(np.ones(3,'d'),grid-1.0)* \
            np.array([max(map(float,coords[0].keys()))-min(map(float,coords[0].keys())),\
                      max(map(float,coords[1].keys()))-min(map(float,coords[1].keys())),\
                      max(map(float,coords[2].keys()))-min(map(float,coords[2].keys())),\
                      ],'d')                                                                        # size from bounding box, corrected for cell-centeredness

  size = np.where(grid > 1, size, min(size[grid > 1]/grid[grid > 1]))                               # spacing for grid==1 equal to smallest among other spacings

# ------------------------------------------ process value field -----------------------------------

  stack = [table.data]

  neighborhood = neighborhoods[options.neighborhood]
  convoluted = np.empty([len(neighborhood)]+list(grid+2),'i')
  microstructure = periodic_3Dpad(np.array(table.data[:,idCol].reshape(grid),'i'))
  
  for i,p in enumerate(neighborhood):
    stencil = np.zeros((3,3,3),'i')
    stencil[1,1,1] = -1
    stencil[p[0]+1,
            p[1]+1,
            p[2]+1] = 1
    convoluted[i,:,:,:] = ndimage.convolve(microstructure,stencil)
  
  distance = np.ones((len(feature_list),grid[0],grid[1],grid[2]),'d')
  
  convoluted = np.sort(convoluted,axis = 0)
  uniques = np.where(convoluted[0,1:-1,1:-1,1:-1] != 0, 1,0)                                        # initialize unique value counter (exclude myself [= 0])

  for i in xrange(1,len(neighborhood)):                                                             # check remaining points in neighborhood
    uniques += np.where(np.logical_and(
                           convoluted[i,1:-1,1:-1,1:-1] != convoluted[i-1,1:-1,1:-1,1:-1],          # flip of ID difference detected?
                           convoluted[i,1:-1,1:-1,1:-1] != 0),                                      # not myself?
                           1,0)                                                                     # count flip

  for i,feature_id in enumerate(feature_list):
    distance[i,:,:,:] = np.where(uniques >= features[feature_id]['aliens'],0.0,1.0)                 # seed with 0.0 when enough unique neighbor IDs are present
    distance[i,:,:,:] = ndimage.morphology.distance_transform_edt(distance[i,:,:,:])*[options.scale]*3

  distance.shape = ([len(feature_list),grid.prod(),1])
  for i in xrange(len(feature_list)):
    stack.append(distance[i,:])

# ------------------------------------------ output result -----------------------------------------

  if len(stack) > 1: table.data = np.hstack(tuple(stack))
  table.data_writeArray('%.12g')

# ------------------------------------------ output finalization -----------------------------------

  table.close()                                                                                     # close input ASCII table (works for stdin)
