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

#--------------------------------------------------------------------------------------------------
#                                MAIN
#--------------------------------------------------------------------------------------------------
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
Produce geom files containing Euclidean distance to grain structural features:
boundaries, triple lines, and quadruple points.

""", version = scriptID)

parser.add_option('-t','--type',          dest = 'type', action = 'extend', type = 'string', metavar = '<string LIST>',
                  help = 'feature type (%s) '%(', '.join(map(lambda x:'|'.join(x['names']),features))) )
parser.add_option('-n','--neighborhood',  dest='neighborhood', choices = neighborhoods.keys(), metavar = 'string',
                  help = 'type of neighborhood (%s) [neumann]'%(', '.join(neighborhoods.keys())))
parser.add_option('-s', '--scale',        dest = 'scale', type = 'float', metavar='float',
                  help = 'voxel size [%default]')

parser.set_defaults(type = [])
parser.set_defaults(neighborhood = 'neumann')
parser.set_defaults(scale = 1.0)

(options,filenames) = parser.parse_args()

if len(options.type) == 0: parser.error('please select a feature type')
if not set(options.type).issubset(set(map(lambda x: x['name'],features))):
  parser.error('type must be chosen from (%s)...'%(', '.join(map(lambda x:', '.join([x['name']]),features))))
if 'biplane' in options.type and 'boundary' in options.type:
  parser.error("please select only one alias for 'biplane' and 'boundary'")
  
feature_list = []
for i,feature in enumerate(features):
  for name in feature['names']:
    for myType in options.type:
      if name.startswith(myType):
        feature_list.append(i)                                                                      # remember valid features
        break

#--- setup file handles ---------------------------------------------------------------------------  
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
                    'output':[open(features[feature]['names'][0]+'_'+name,'w') for feature in feature_list],
                    'croak':sys.stdout,
                    })

#--- loop over input files ------------------------------------------------------------------------
for file in files:
  file['croak'].write('\033[1m' + scriptName + '\033[0m: ' + (file['name'] if file['name'] != 'STDIN' else '') + '\n')

  table = damask.ASCIItable(file['input'],file['output'][0],labels = False)
  table.head_read()

#--- interpret header ----------------------------------------------------------------------------
  info = {
          'grid':    np.zeros(3,'i'),
          'size':    np.zeros(3,'d'),
          'origin':  np.zeros(3,'d'),
          'homogenization':  0,
          'microstructures': 0,
         }
  newInfo = {
          'grid':    np.zeros(3,'i'),
          'origin':  np.zeros(3,'d'),
          'microstructures': 0,
         }
  extra_header = []

  for header in table.info:
    headitems = map(str.lower,header.split())
    if len(headitems) == 0: continue                                                              # skip blank lines
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
  microstructure = np.zeros(info['grid'].prod(),'i')                                            # initialize as flat array
  i = 0

  while table.data_read():
    items = table.data
    if len(items) > 2:
      if   items[1].lower() == 'of': items = [int(items[2])]*int(items[0])
      elif items[1].lower() == 'to': items = xrange(int(items[0]),1+int(items[2]))
      else:                            items = map(int,items)
    else:                              items = map(int,items)

    s = len(items)
    microstructure[i:i+s] = items
    i += s

  
  neighborhood = neighborhoods[options.neighborhood]
  convoluted = np.empty([len(neighborhood)]+list(info['grid']+2),'i')
  structure = periodic_3Dpad(microstructure.reshape(info['grid'],order='F'))
  
  for i,p in enumerate(neighborhood):
    stencil = np.zeros((3,3,3),'i')
    stencil[1,1,1] = -1
    stencil[p[0]+1,
            p[1]+1,
            p[2]+1] = 1
    convoluted[i,:,:,:] = ndimage.convolve(structure,stencil)
  
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

  for i,feature in enumerate(feature_list):
    newInfo['microstructures'] = int(math.ceil(distance[i,:,:,:].max()))

#--- write header ---------------------------------------------------------------------------------
    table = damask.ASCIItable(file['input'],file['output'][i],labels = False)
    table.labels_clear()
    table.info_clear()
    table.info_append(extra_header+[
      scriptID + ' ' + ' '.join(sys.argv[1:]),
      "grid\ta %i\tb %i\tc %i"%(info['grid'][0],info['grid'][1],info['grid'][2],),
      "size\tx %f\ty %f\tz %f"%(info['size'][0],info['size'][1],info['size'][2],),
      "origin\tx %f\ty %f\tz %f"%(info['origin'][0],info['origin'][1],info['origin'][2],),
      "homogenization\t%i"%info['homogenization'],
      "microstructures\t%i"%(newInfo['microstructures']),
      ])
    table.head_write()
    table.output_flush()
    
# --- write microstructure information ------------------------------------------------------------
    formatwidth = int(math.floor(math.log10(distance[i,:,:,:].max())+1))
    table.data = distance[i,:,:,:].reshape((info['grid'][0],info['grid'][1]*info['grid'][2]),order='F').transpose()
    table.data_writeArray('%%%ii'%(formatwidth),delimiter=' ')
    file['output'][i].close()
    
#--- output finalization --------------------------------------------------------------------------
  if file['name'] != 'STDIN':
    table.input_close()  
