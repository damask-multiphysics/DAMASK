#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,re,sys,math,string
import numpy as np
from optparse import OptionParser
from scipy import ndimage
import damask

scriptID   = string.replace('$Id$','\n','\\n')
scriptName = scriptID.split()[1]

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

features = [ \
            {'aliens': 1, 'name': 'biplane'},
            {'aliens': 1, 'name': 'boundary'},
            {'aliens': 2, 'name': 'tripleline'},
            {'aliens': 3, 'name': 'quadruplepoint'}
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
                                        [-1,-1, 0],
                                        [ 0,-1, 0],
                                        [ 1,-1, 0],
                                        [-1, 0, 0],
#
                                        [ 1, 0, 0],
                                        [-1, 1, 0],
                                        [ 0, 1, 0],
                                        [ 1, 1, 0],
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

parser.add_option('-c','--coordinates', dest='coords', action='store', type='string', metavar='string',
                                        help='column heading for coordinates [%default]')
parser.add_option('-i','--identifier',  dest='id', action='store', type='string', metavar = 'string',
                                        help='heading of column containing grain identifier [%default]')
parser.add_option('-t','--type',        dest='type', action='extend', type='string', metavar='<string LIST>',
                                        help='feature type (%s)'%(', '.join(map(lambda x:', '.join([x['name']]),features))))
parser.add_option('-n','--neighborhood',dest='neigborhood', action='store', type='choice', 
                                        choices=neighborhoods.keys(), metavar='string',
                                        help='type of neighborhood (%s) [neumann]'%(', '.join(neighborhoods.keys())))
parser.set_defaults(type = [])
parser.set_defaults(coords = 'ip')
parser.set_defaults(id = 'texture')
parser.set_defaults(neighborhood = 'neumann')

(options,filenames) = parser.parse_args()

if len(options.type) == 0: parser.error('please select a feature type')
if not set(options.type).issubset(set(map(lambda x: x['name'],features))):
  parser.error('type must be chosen from (%s)...'%(', '.join(map(lambda x:', '.join([x['name']]),features))))
if 'biplane' in options.type and 'boundary' in options.type:
  parser.error("please select only one alias for 'biplane' and 'boundary'")

feature_list = []
for i,feature in enumerate(features):
  if feature['name'] in options.type: feature_list.append(i)                                        # remember valid features
# ------------------------------------------ setup file handles -----------------------------------

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
    table.labels_append('ED_%s(%s)'%(features[feature]['name'],options.id))                         # extend ASCII header with new labels

  table.head_write()

# ------------------------------------------ process data ---------------------------------------   
  
  table.data_readArray([options.coords+'.x',options.coords+'.y',options.coords+'.z',options.id])

  grid = [{},{},{}]
  for i in xrange(len(table.data)):
    for j in xrange(3):
      grid[j][str(table.data[i,j])] = True

  resolution = np.array(map(len,grid),'i')
  unitlength = 0.0
  for i,r in enumerate(resolution):
    if r > 1: unitlength = max(unitlength,(max(map(float,grid[i].keys()))-min(map(float,grid[i].keys())))/(r-1.0))

  neighborhood = neighborhoods[options.neighborhood]
  convoluted = np.empty([len(neighborhood)]+list(resolution+2),'i')
  microstructure = periodic_3Dpad(np.array(table.data[:,3].reshape(resolution),'i'))
  
  for i,p in enumerate(neighborhood):
    stencil = np.zeros((3,3,3),'i')
    stencil[1,1,1] = -1
    stencil[p[0]+1,
            p[1]+1,
            p[2]+1] = 1

    convoluted[i,:,:,:] = ndimage.convolve(microstructure,stencil)
  
  distance = np.ones((len(feature_list),resolution[0],resolution[1],resolution[2]),'d')
  
  convoluted = np.sort(convoluted,axis=0)
  uniques = np.zeros(resolution)
  check = np.empty(resolution)
  check[:,:,:] = np.nan
  for i in xrange(len(neighborhood)):
    uniques += np.where(convoluted[i,1:-1,1:-1,1:-1] == check,0,1)
    check = convoluted[i,1:-1,1:-1,1:-1]
  for i,feature_id in enumerate(feature_list):
    distance[i,:,:,:] = np.where(uniques > features[feature_id]['aliens'],0.0,1.0)
  
  for i in xrange(len(feature_list)):
    distance[i,:,:,:] = ndimage.morphology.distance_transform_edt(distance[i,:,:,:])*[unitlength]*3
  distance.shape = (len(feature_list),resolution.prod())
  
# ------------------------------------------ process data ---------------------------------------
  table.data_rewind()
  l = 0
  while table.data_read():
    for i in xrange(len(feature_list)):
      table.data_append(distance[i,l])                                                              # add all distance fields
    l += 1
    outputAlive = table.data_write()                                                                # output processed line

# ------------------------------------------ output result ---------------------------------------  
  outputAlive and table.output_flush()                                                              # just in case of buffered ASCII table

  file['input'].close()                                                                             # close input ASCII table (works for stdin)
  file['output'].close()                                                                            # close output ASCII table (works for stdout)
  os.rename(file['name']+'_tmp',file['name'])                                                       # overwrite old one with tmp new
