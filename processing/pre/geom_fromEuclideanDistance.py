#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,re,sys,math,numpy,string,damask
from scipy import ndimage
from optparse import OptionParser, OptionGroup, Option, SUPPRESS_HELP

scriptID = '$Id$'
scriptName = scriptID.split()[1]

#--------------------------------------------------------------------------------------------------
class extendableOption(Option):
#--------------------------------------------------------------------------------------------------
# used for definition of new option parser action 'extend', which enables to take multiple option arguments
# taken from online tutorial http://docs.python.org/library/optparse.html
  
  ACTIONS = Option.ACTIONS + ("extend",)
  STORE_ACTIONS = Option.STORE_ACTIONS + ("extend",)
  TYPED_ACTIONS = Option.TYPED_ACTIONS + ("extend",)
  ALWAYS_TYPED_ACTIONS = Option.ALWAYS_TYPED_ACTIONS + ("extend",)

  def take_action(self, action, dest, opt, value, values, parser):
    if action == "extend":
      lvalue = value.split(",")
      values.ensure_value(dest, []).extend(lvalue)
    else:
      Option.take_action(self, action, dest, opt, value, values, parser)


def periodic_3Dpad(array, rimdim=(1,1,1)):

  rimdim = numpy.array(rimdim,'i')
  size = numpy.array(array.shape,'i')
  padded = numpy.empty(size+2*rimdim,array.dtype)
  padded[rimdim[0]:rimdim[0]+size[0],
         rimdim[1]:rimdim[1]+size[1],
         rimdim[2]:rimdim[2]+size[2]] = array

  p = numpy.zeros(3,'i')
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

features = [
            {'aliens': 1, 'names': ['boundary','biplane'],},
            {'aliens': 2, 'names': ['tripleline',],},
            {'aliens': 3, 'names': ['quadruplepoint',],}
           ]

neighborhoods = {
                  'neumann':numpy.array([
                                        [-1, 0, 0],
                                        [ 1, 0, 0],
                                        [ 0,-1, 0],
                                        [ 0, 1, 0],
                                        [ 0, 0,-1],
                                        [ 0, 0, 1],
                                      ]),
                  'moore':numpy.array([
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

parser = OptionParser(option_class=extendableOption, usage='%prog options [file[s]]', description = """
Produce geom files containing Euclidean distance to grain structural features:
boundaries, triple lines, and quadruple points.
""" + string.replace(scriptID,'\n','\\n')
)

parser.add_option('-t','--type',          dest = 'type', action = 'extend', type = 'string', metavar = '<string LIST>',
                  help = 'feature type (%s) '%(', '.join(map(lambda x:'|'.join(x['names']),features))) )
parser.add_option('-n','--neighborhood',  dest='neighborhood', choices = neighborhoods.keys(), metavar = 'string',
                  help = 'type of neighborhood (%s) [neumann]'%(', '.join(neighborhoods.keys())))
parser.add_option('-s', '--scale',        dest = 'scale', type = 'float',
                  help = 'voxel size [%default]')

parser.set_defaults(type = [])
parser.set_defaults(neighborhood = 'neumann')
parser.set_defaults(scale = 1.0)

(options,filenames) = parser.parse_args()

  
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
  if file['name'] != 'STDIN': file['croak'].write('\033[1m'+scriptName+'\033[0m: '+file['name']+'\n')
  else: file['croak'].write('\033[1m'+scriptName+'\033[0m\n')

  theTable = damask.ASCIItable(file['input'],file['output'][0],labels = False)
  theTable.head_read()

#--- interpret header ----------------------------------------------------------------------------
  info = {
          'grid':    numpy.zeros(3,'i'),
          'size':    numpy.zeros(3,'d'),
          'origin':  numpy.zeros(3,'d'),
          'homogenization':  0,
          'microstructures': 0,
         }
  newInfo = {
          'grid':    numpy.zeros(3,'i'),
          'origin':  numpy.zeros(3,'d'),
          'microstructures': 0,
         }
  extra_header = []

  for header in theTable.info:
    headitems = map(str.lower,header.split())
    if len(headitems) == 0: continue                                                              # skip blank lines
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

  if numpy.any(info['grid'] < 1):
    file['croak'].write('invalid grid a b c.\n')
    continue
  if numpy.any(info['size'] <= 0.0):
    file['croak'].write('invalid size x y z.\n')
    continue

#--- read data ------------------------------------------------------------------------------------
  microstructure = numpy.zeros(info['grid'].prod(),'i')                                            # initialize as flat array
  i = 0

  while theTable.data_read():
    items = theTable.data
    if len(items) > 2:
      if   items[1].lower() == 'of': items = [int(items[2])]*int(items[0])
      elif items[1].lower() == 'to': items = xrange(int(items[0]),1+int(items[2]))
      else:                            items = map(int,items)
    else:                              items = map(int,items)

    s = len(items)
    microstructure[i:i+s] = items
    i += s

  
  neighborhood = neighborhoods[options.neighborhood]
  convoluted = numpy.empty([len(neighborhood)]+list(info['grid']+2),'i')
  structure = periodic_3Dpad(microstructure.reshape(info['grid'],order='F'))
  
  for i,p in enumerate(neighborhood):
    stencil = numpy.zeros((3,3,3),'i')
    stencil[1,1,1] = -1
    stencil[p[0]+1,
            p[1]+1,
            p[2]+1] = 1
    convoluted[i,:,:,:] = ndimage.convolve(structure,stencil)
  
  distance = numpy.ones((len(feature_list),info['grid'][0],info['grid'][1],info['grid'][2]),'d')
  
  convoluted = numpy.sort(convoluted,axis = 0)
  uniques = numpy.where(convoluted[0,1:-1,1:-1,1:-1] != 0, 1,0)                   # initialize unique value counter (exclude myself [= 0])

  for i in xrange(1,len(neighborhood)):                                           # check remaining points in neighborhood
    uniques += numpy.where(numpy.logical_and(
                           convoluted[i,1:-1,1:-1,1:-1] != convoluted[i-1,1:-1,1:-1,1:-1],    # flip of ID difference detected?
                           convoluted[i,1:-1,1:-1,1:-1] != 0),                                # not myself?
                           1,0)                                                   # count flip

  for i,feature_id in enumerate(feature_list):
    distance[i,:,:,:] = numpy.where(uniques >= features[feature_id]['aliens'],0.0,1.0) # seed with 0.0 when enough unique neighbor IDs are present

  for i in xrange(len(feature_list)):
    distance[i,:,:,:] = ndimage.morphology.distance_transform_edt(distance[i,:,:,:])*[options.scale]*3

  for i,feature in enumerate(feature_list):
    newInfo['microstructures'] = int(math.ceil(distance[i,:,:,:].max()))

#--- write header ---------------------------------------------------------------------------------
    theTable = damask.ASCIItable(file['input'],file['output'][i],labels = False)
    theTable.labels_clear()
    theTable.info_clear()
    theTable.info_append(extra_header+[
      scriptID + ' ' + ' '.join(sys.argv[1:]),
      "grid\ta %i\tb %i\tc %i"%(info['grid'][0],info['grid'][1],info['grid'][2],),
      "size\tx %f\ty %f\tz %f"%(info['size'][0],info['size'][1],info['size'][2],),
      "origin\tx %f\ty %f\tz %f"%(info['origin'][0],info['origin'][1],info['origin'][2],),
      "homogenization\t%i"%info['homogenization'],
      "microstructures\t%i"%(newInfo['microstructures']),
      ])
    theTable.head_write()
    theTable.output_flush()
    
# --- write microstructure information ------------------------------------------------------------
    formatwidth = int(math.floor(math.log10(distance[i,:,:,:].max())+1))
    theTable.data = distance[i,:,:,:].reshape((info['grid'][0],info['grid'][1]*info['grid'][2]),order='F').transpose()
    theTable.data_writeArray('%%%ii'%(formatwidth),delimiter=' ')
    file['output'][i].close()
    
#--- output finalization --------------------------------------------------------------------------
  if file['name'] != 'STDIN':
    file['input'].close()
