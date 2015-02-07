#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,string,re,math
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
parser.add_option('-t','--type',          dest = 'type', action = 'extend', type = 'string', metavar = '<string LIST>',
                  help = 'feature type (%s) '%(', '.join(map(lambda x:'|'.join(x['names']),features))) )
parser.add_option('-n','--neighborhood',  dest='neighborhood', choices = neighborhoods.keys(), metavar = 'string',
                  help = 'type of neighborhood (%s) [neumann]'%(', '.join(neighborhoods.keys())))
parser.add_option('-s', '--scale',        dest = 'scale', type = 'float',
                  help = 'voxel size [%default]')
parser.set_defaults(type = [])
parser.set_defaults(coords = 'ip')
parser.set_defaults(id = 'texture')
parser.set_defaults(neighborhood = 'neumann')
parser.set_defaults(scale = 1.0)

(options,filenames) = parser.parse_args()

if len(options.type) == 0: parser.error('please select a feature type')
if 'biplane' in options.type and 'boundary' in options.type:
  parser.error("please select only one alias for 'biplane' and 'boundary'")

feature_list = []
for i,feature in enumerate(features):
  for name in feature['names']:
    for myType in options.type:
      if name.startswith(myType):
        feature_list.append(i)                                                                      # remember valid features
        break

files = []
for name in filenames:
  if os.path.exists(name):
    files.append({'name':name, 'input':open(name), 'output':open(name+'_tmp','w'), 'croak':sys.stderr})

# ------------------------------------------ loop over input files ---------------------------------
for file in files:
  file['croak'].write('\033[1m'+scriptName+'\033[0m: '+file['name']+'\n')

  table = damask.ASCIItable(file['input'],file['output'],False)                                     # make unbuffered ASCII_table
  table.head_read()                                                                                 # read ASCII header info
  table.info_append(scriptID + '\t' + ' '.join(sys.argv[1:]))

# --------------- figure out position of labels and coordinates ------------------------------------
  try:
    locationCol = table.labels.index('%s.x'%options.coords)                                         # columns containing location data
  except ValueError:
    file['croak'].write('no coordinate data (%s.x) found...\n'%options.coords)
    continue

  if options.id not in table.labels:
    file['croak'].write('column %s not found...\n'%options.id)
    continue

# ------------------------------------------ assemble header ---------------------------------------
  for feature in feature_list:
    table.labels_append('ED_%s(%s)'%(features[feature]['names'],options.id))                         # extend ASCII header with new labels

  table.head_write()

# ------------------------------------------ process data ------------------------------------------
  
  table.data_readArray(['1_'+options.coords,'2_'+options.coords,'3_'+options.coords,options.id])

  coords = [{},{},{}]
  for i in xrange(len(table.data)):
    for j in xrange(3):
      coords[j][str(table.data[i,j])] = True

  grid = np.array(map(len,coords),'i')
  unitlength = 0.0
  for i,r in enumerate(grid):
    if r > 1: unitlength = max(unitlength,(max(map(float,coords[i].keys()))-min(map(float,coords[i].keys())))/(r-1.0))

  neighborhood = neighborhoods[options.neighborhood]
  convoluted = np.empty([len(neighborhood)]+list(grid+2),'i')
  microstructure = periodic_3Dpad(np.array(table.data[:,3].reshape(grid),'i'))
  
  for i,p in enumerate(neighborhood):
    stencil = np.zeros((3,3,3),'i')
    stencil[1,1,1] = -1
    stencil[p[0]+1,
            p[1]+1,
            p[2]+1] = 1
    convoluted[i,:,:,:] = ndimage.convolve(microstructure,stencil)
  
  distance = np.ones((len(feature_list),info['grid'][0],info['grid'][1],info['grid'][2]),'d')
  
  convoluted = np.sort(convoluted,axis = 0)
  uniques = np.where(convoluted[0,1:-1,1:-1,1:-1] != 0, 1,0)                   # initialize unique value counter (exclude myself [= 0])

  for i in xrange(1,len(neighborhood)):                                           # check remaining points in neighborhood
    uniques += np.where(np.logical_and(
                           convoluted[i,1:-1,1:-1,1:-1] != convoluted[i-1,1:-1,1:-1,1:-1],    # flip of ID difference detected?
                           convoluted[i,1:-1,1:-1,1:-1] != 0),                                # not myself?
                           1,0)                                                   # count flip

  for i,feature_id in enumerate(feature_list):
    distance[i,:,:,:] = np.where(uniques >= features[feature_id]['aliens'],0.0,1.0) # seed with 0.0 when enough unique neighbor IDs are present

  for i in xrange(len(feature_list)):
    distance[i,:,:,:] = ndimage.morphology.distance_transform_edt(distance[i,:,:,:])*[options.scale]*3
  distance.shape = (len(feature_list),grid.prod())
  
# ------------------------------------------ process data ------------------------------------------
  table.data_rewind()
  l = 0
  while table.data_read():
    for i in xrange(len(feature_list)):
      table.data_append(distance[i,l])                                                              # add all distance fields
    l += 1
    outputAlive = table.data_write()                                                                # output processed line

# ------------------------------------------ output result -----------------------------------------
  outputAlive and table.output_flush()                                                              # just in case of buffered ASCII table

  table.input_close()                                                                               # close input ASCII table
  table.output_close()                                                                              # close output ASCII table
  os.rename(file['name']+'_tmp',file['name'])                                                       # overwrite old one with tmp new
