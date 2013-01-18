#!/usr/bin/env python

import os,re,sys,math,numpy,string,damask
from scipy import ndimage
from optparse import OptionParser, Option

# -----------------------------
class extendableOption(Option):
# -----------------------------
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

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

features = [ \
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
Add column(s) containing Euclidean distance to grain structural features:
boundaries, triple lines, and quadruple points.

""" + string.replace('$Id$','\n','\\n')
)

parser.add_option('-i','--identifier',  dest='id', action='store', type='string', \
                                        help='heading of column containing grain identifier [%default]', \
                                        metavar='<label>')
parser.add_option('-t','--type',        dest='type', action='extend', type='string', \
                                        help='feature type (%s)'%(', '.join(map(lambda x:', '.join(x['names']),features))))
parser.add_option('-n','--neighborhood',  dest='neigborhood', action='store', type='string', \
                                        help='type of neighborhood (%s)'%(', '.join(neighborhoods.keys())), \
                                        metavar='<int>')

parser.set_defaults(type = [])
parser.set_defaults(id = 'texture')
parser.set_defaults(neighborhood = 'neumann')

(options,filenames) = parser.parse_args()

options.neighborhood = options.neighborhood.lower()
if options.neighborhood not in neighborhoods:
  parser.error('unknown neighborhood %s!'%options.neighborhood)
  
feature_list = []
for i,feature in enumerate(features):
  for name in feature['names']:
    for type in options.type:
      if name.startswith(type):
        feature_list.append(i)                            # remember valid features
        break

# ------------------------------------------ setup file handles ---------------------------------------  

files = []
if filenames == []:
  files.append({'name':'STDIN', 'input':sys.stdin, 'output':sys.stdout, 'croak':sys.stderr})
else:
  for name in filenames:
    if os.path.exists(name):
      files.append({'name':name, 'input':open(name), 'output':open(name+'_tmp','w'), 'croak':sys.stderr})

# ------------------------------------------ loop over input files ---------------------------------------  

for file in files:
  if file['name'] != 'STDIN': file['croak'].write(file['name']+'\n')

  table = damask.ASCIItable(file['input'],file['output'],False)             # make unbuffered ASCII_table
  table.head_read()                                                         # read ASCII header info
  table.info_append(string.replace('$Id$','\n','\\n') + \
                    '\t' + ' '.join(sys.argv[1:]))

# ------------------------------------------ assemble header ---------------------------------------  

  if options.id not in table.labels:
    file['croak'].write('column %s not found...\n'%options.id)
    continue

  for feature in feature_list:
    table.labels_append('ED_%s(%s)'%(features[feature]['names'][0],options.id))     # extend ASCII header with new labels

  table.head_write()

# ------------------------------------------ process data ---------------------------------------   
  
  structure = table.data_asArray(['ip.x','ip.y','ip.z',options.id])

  grid = [{},{},{}]
  for i in xrange(len(structure)):
    for j in xrange(3):
      grid[j][str(structure[i,j])] = True

  resolution = numpy.array(map(len,grid),'i')
  unitlength = 0.0
  for i,r in enumerate(resolution):
    if r > 1: unitlength = max(unitlength,(max(map(float,grid[i].keys()))-min(map(float,grid[i].keys())))/(r-1.0))

  neighborhood = neighborhoods[options.neighborhood]
  convoluted = numpy.empty([len(neighborhood)]+list(resolution+2),'i')
  microstructure = periodic_3Dpad(numpy.array(structure[:,3].reshape(resolution),'i'))
  
  for i,p in enumerate(neighborhood):
    stencil = numpy.zeros((3,3,3),'i')
    stencil[1,1,1] = -1
    stencil[p[0]+1,
            p[1]+1,
            p[2]+1] = 1

    convoluted[i,:,:,:] = ndimage.convolve(microstructure,stencil)
  
  distance = numpy.ones((len(feature_list),resolution[0],resolution[1],resolution[2]),'d')
  
  convoluted = numpy.sort(convoluted,axis=0)
  uniques = numpy.zeros(resolution)
  check = numpy.empty(resolution)
  check[:,:,:] = numpy.nan
  for i in xrange(len(neighborhood)):
    uniques += numpy.where(convoluted[i,1:-1,1:-1,1:-1] == check,0,1)
    check = convoluted[i,1:-1,1:-1,1:-1]
  for i,feature_id in enumerate(feature_list):
    distance[i,:,:,:] = numpy.where(uniques > features[feature_id]['aliens'],0.0,1.0)
  
  for i in xrange(len(feature_list)):
    distance[i,:,:,:] = ndimage.morphology.distance_transform_edt(distance[i,:,:,:])*[unitlength]*3
  distance.shape = (len(feature_list),resolution.prod())
  
  table.data_rewind()
  l = 0
  while table.data_read():
    for i in xrange(len(feature_list)):
      table.data_append(distance[i,l])                                      # add all distance fields
    table.data_write()                                                      # output processed line
    l += 1
    
# ------------------------------------------ output result ---------------------------------------  

  table.output_flush()                                                      # just in case of buffered ASCII table

  if file['name'] != 'STDIN':
    file['input'].close()                                                   # close input ASCII table
    file['output'].close()                                                  # close output ASCII table
    os.rename(file['name']+'_tmp',file['name'])                             # overwrite old one with tmp new
    