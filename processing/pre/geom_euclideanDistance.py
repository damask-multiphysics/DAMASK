#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,re,sys,math,numpy,string,damask
from scipy import ndimage
from optparse import OptionParser, Option

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
            {'aliens': 1, 'names': ['boundary, biplane'],},
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

parser = OptionParser(option_class=extendableOption, usage='%prog options [file[s]]', description = """
Produce geom files containing Euclidean distance to grain structural features:
boundaries, triple lines, and quadruple points.

""" + string.replace('$Id$','\n','\\n')
)

parser.add_option('-t','--type',        dest='type', action='extend', type='string', \
                  help='feature type (%s)'%(', '.join(map(lambda x:', '.join(x['names']),features))))
parser.add_option('-n','--neighborhood',  dest='neigborhood', action='store', type='string', \
                  help='type of neighborhood (%s) [neumann]'%(', '.join(neighborhoods.keys())))
parser.add_option('-2', '--twodimensional', dest='twoD', action='store_true', \
                  help='output geom file with two-dimensional data arrangement [%default]')

parser.set_defaults(type = [])
parser.set_defaults(neighborhood = 'neumann')
parser.set_defaults(twoD = False)

(options,filenames) = parser.parse_args()

options.neighborhood = options.neighborhood.lower()
if options.neighborhood not in neighborhoods:
  parser.error('unknown neighborhood %s!'%options.neighborhood)
  
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
  if file['name'] != 'STDIN': file['croak'].write(file['name']+'\n')

  firstline = file['input'].readline()
  m = re.search('(\d+)\s*head', firstline.lower())
  if m:
    headerlines = int(m.group(1))
    headers  = [file['input'].readline() for i in range(headerlines)]
  else:
    headerlines = 1
    headers = firstline

  content = file['input'].readlines()
  file['input'].close()

#--- interpretate header --------------------------------------------------------------------------
  info = {
          'grid':   numpy.array([0,0,0]),
          'size':   numpy.array([0.0,0.0,0.0]),
          'origin': numpy.array([0.0,0.0,0.0]),
          'microstructures': 0,
          'homogenization':  0
         }

  newInfo = {
        	'microstructures': 0,
            }

  new_header = []
  new_header.append('$Id$\n')
  for header in headers:
    headitems = map(str.lower,header.split())
    if headitems[0] == 'resolution': headitems[0] = 'grid'
    if headitems[0] == 'dimension':  headitems[0] = 'size'
    if headitems[0] in mappings.keys():
      if headitems[0] in identifiers.keys():
        for i in xrange(len(identifiers[headitems[0]])):
          info[headitems[0]][i] = \
            mappings[headitems[0]](headitems[headitems.index(identifiers[headitems[0]][i])+1])
      else:
        info[headitems[0]] = mappings[headitems[0]](headitems[1])
    else:
      new_header.append(header)

  if numpy.all(info['grid'] == 0):
    file['croak'].write('no grid info found.\n')
    continue
  if numpy.all(info['size'] == 0.0):
    file['croak'].write('no size info found.\n')
    continue

  file['croak'].write('grid     a b c:  %s\n'%(' x '.join(map(str,info['grid']))) + \
                      'size     x y z:  %s\n'%(' x '.join(map(str,info['size']))) + \
                      'origin   x y z:  %s\n'%(' : '.join(map(str,info['origin']))) + \
                      'homogenization:  %i\n'%info['homogenization'] + \
                      'microstructures: %i\n'%info['microstructures'])

  new_header.append("grid\ta %i\tb %i\tc %i\n"%(info['grid'][0],info['grid'][1],info['grid'][2],))
  new_header.append("size\tx %f\ty %f\tz %f\n"%(info['size'][0],info['size'][1],info['size'][2],))
  new_header.append("origin\tx %f\ty %f\tz %f\n"%(info['origin'][0],info['origin'][1],info['origin'][2],))
  new_header.append("homogenization\t%i\n"%info['homogenization'])

#--- process input --------------------------------------------------------------------------------
  structure = numpy.zeros(info['grid'],'i')
  i = 0
  for line in content:  
    for item in map(int,line.split()):
      structure[i%info['grid'][0],
               (i/info['grid'][0])%info['grid'][1],
                i/info['grid'][0] /info['grid'][1]] = item
      i += 1
  
  neighborhood = neighborhoods[options.neighborhood]
  convoluted = numpy.empty([len(neighborhood)]+list(info['grid']+2),'i')
  microstructure = periodic_3Dpad(structure)
  
  for i,p in enumerate(neighborhood):
    stencil = numpy.zeros((3,3,3),'i')
    stencil[1,1,1] = -1
    stencil[p[0]+1,
            p[1]+1,
            p[2]+1] = 1

    convoluted[i,:,:,:] = ndimage.convolve(microstructure,stencil)
  
  distance = numpy.ones((len(feature_list),info['grid'][0],info['grid'][1],info['grid'][2]),'d')
  
  convoluted = numpy.sort(convoluted,axis=0)
  uniques = numpy.zeros(info['grid'])
  check = numpy.empty(info['grid'])
  check[:,:,:] = numpy.nan
  for i in xrange(len(neighborhood)):
    uniques += numpy.where(convoluted[i,1:-1,1:-1,1:-1] == check,0,1)
    check = convoluted[i,1:-1,1:-1,1:-1]
  for i,feature_id in enumerate(feature_list):
    distance[i,:,:,:] = numpy.where(uniques > features[feature_id]['aliens'],0.0,1.0)
  
  for i in xrange(len(feature_list)):
    distance[i,:,:,:] = ndimage.morphology.distance_transform_edt(distance[i,:,:,:])*\
                                                [max(info['size']/info['grid'])]*3
  for i,feature in enumerate(feature_list):
    newInfo['microstructures'] = int(math.ceil(distance[i,:,:,:].max()))
    formatwidth = int(math.floor(math.log10(distance[i,:,:,:].max())+1))

#--- assemble header and report changes -----------------------------------------------------------
    output  = '%i\theader\n'%(len(new_header)+1)
    output += ''.join(new_header)
    output += "microstructures\t%i\n"%newInfo['microstructures']
    file['croak'].write('\n'+features[i]['names'][0]+'\n')
    if (newInfo['microstructures'] != info['microstructures']):
      file['croak'].write('--> microstructures: %i\n'%newInfo['microstructures'])

#--- write new data -------------------------------------------------------------------------------
    for z in xrange(info['grid'][2]):
      for y in xrange(info['grid'][1]):
        output += {True:' ',False:'\n'}[options.twoD].join(map(lambda x: \
                                       ('%%%ii'%formatwidth)%(round(x)), distance[i,:,y,z])) + '\n'
    file['output'][i].write(output)
	
    if file['name'] != 'STDIN':
      file['output'][i].close()

#--- output finalization -------------------------------------------------------------------------- 
  if file['name'] != 'STDIN':
    file['input'].close()
