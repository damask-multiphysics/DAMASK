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

identifiers = {
        'resolution': ['a','b','c'],
        'dimension':  ['x','y','z'],
        'origin':     ['x','y','z'],
          }
mappings = {
        'resolution': lambda x: int(x),
        'dimension':  lambda x: float(x),
        'origin':     lambda x: float(x),
          }

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
Produce geom files containing Euclidean distance to grain structural features:
boundaries, triple lines, and quadruple points.

""" + string.replace('$Id: addEuclideanDistance.py 2039 2012-12-19 14:50:45Z MPIE\p.shanthraj $','\n','\\n')
)

parser.add_option('-t','--type',        dest='type', action='extend', type='string', \
                                        help='feature type (%s)'%(', '.join(map(lambda x:', '.join(x['names']),features))))
parser.add_option('-n','--neighborhood',  dest='neigborhood', action='store', type='string', \
                                        help='type of neighborhood (%s)'%(', '.join(neighborhoods.keys())), \
                                        metavar='<int>')
parser.add_option('-2', '--twodimensional', dest='twoD', action='store_true', \
                  help='output geom file with two-dimensional data arrangement')

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
    for type in options.type:
      if name.startswith(type):
        feature_list.append(i)                            # remember valid features
        break

print feature_list
# ------------------------------------------ setup file handles ---------------------------------------  

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

# ------------------------------------------ loop over input files ---------------------------------------  

for file in files:
  if file['name'] != 'STDIN': file['croak'].write(file['name']+'\n')

  #  get labels by either read the first row, or - if keyword header is present - the last line of the header

  firstline = file['input'].readline()
  m = re.search('(\d+)\s*head', firstline.lower())
  if m:
    headerlines = int(m.group(1))
    headers  = [firstline]+[file['input'].readline() for i in range(headerlines)]
  else:
    headerlines = 1
    headers = firstline

  content = file['input'].readlines()
  file['input'].close()

  info = {'resolution': numpy.array([0,0,0]),
          'dimension':  numpy.array([0.0,0.0,0.0]),
          'origin':     numpy.array([0.0,0.0,0.0]),
          'homogenization': 1,
         }

  new_header = []
  for header in headers:
    headitems = map(str.lower,header.split())
    if headitems[0] in mappings.keys():
      if headitems[0] in identifiers.keys():
        for i in xrange(len(identifiers[headitems[0]])):
          info[headitems[0]][i] = \
            mappings[headitems[0]](headitems[headitems.index(identifiers[headitems[0]][i])+1])
      else:
        info[headitems[0]] = mappings[headitems[0]](headitems[1])

  if numpy.all(info['resolution'] == 0):
    file['croak'].write('no resolution info found.\n')
    continue
  if numpy.all(info['dimension'] == 0.0):
    file['croak'].write('no dimension info found.\n')
    continue

  file['croak'].write('resolution:     %s\n'%(' x '.join(map(str,info['resolution']))) + \
                      'dimension:      %s\n'%(' x '.join(map(str,info['dimension']))) + \
                      'origin:         %s\n'%(' : '.join(map(str,info['origin']))) + \
                      'homogenization: %i\n'%info['homogenization'])

  new_header.append("resolution\ta %i\tb %i\tc %i\n"%( 
    info['resolution'][0],
    info['resolution'][1],
    info['resolution'][2],))
  new_header.append("dimension\tx %f\ty %f\tz %f\n"%(
    info['dimension'][0],
    info['dimension'][1],
    info['dimension'][2],))
  new_header.append("origin\tx %f\ty %f\tz %f\n"%(
    info['origin'][0],
    info['origin'][1],
    info['origin'][2],))
  new_header.append("homogenization\t%i\n"%info['homogenization'])

  structure = numpy.zeros(info['resolution'],'i')
  i = 0
  for line in content:  
    for item in map(int,line.split()):
      structure[i%info['resolution'][0],
                (i/info['resolution'][0])%info['resolution'][1],
                 i/info['resolution'][0] /info['resolution'][1]] = item
      i += 1
  
  neighborhood = neighborhoods[options.neighborhood]
  convoluted = numpy.empty([len(neighborhood)]+list(info['resolution']+2),'i')
  microstructure = periodic_3Dpad(structure)
  
  for i,p in enumerate(neighborhood):
    stencil = numpy.zeros((3,3,3),'i')
    stencil[1,1,1] = -1
    stencil[p[0]+1,
            p[1]+1,
            p[2]+1] = 1

    convoluted[i,:,:,:] = ndimage.convolve(microstructure,stencil)
  
  distance = numpy.ones((len(feature_list),info['resolution'][0],info['resolution'][1],info['resolution'][2]),'d')
  
  convoluted = numpy.sort(convoluted,axis=0)
  uniques = numpy.zeros(info['resolution'])
  check = numpy.empty(info['resolution'])
  check[:,:,:] = numpy.nan
  for i in xrange(len(neighborhood)):
    uniques += numpy.where(convoluted[i,1:-1,1:-1,1:-1] == check,0,1)
    check = convoluted[i,1:-1,1:-1,1:-1]
  for i,feature_id in enumerate(feature_list):
    distance[i,:,:,:] = numpy.where(uniques > features[feature_id]['aliens'],0.0,1.0)
  
  for i in xrange(len(feature_list)):
    distance[i,:,:,:] = ndimage.morphology.distance_transform_edt(distance[i,:,:,:])*[max(info['dimension']/info['resolution'])]*3


  for i,feature in enumerate(feature_list):
		formatwidth = int(math.floor(math.log10(distance[i,:,:,:].max())+1))
            
# ------------------------------------------ assemble header ---------------------------------------  

		output  = '%i\theader\n'%(len(new_header))
		output += ''.join(new_header)

# ------------------------------------- regenerate texture information ----------------------------------  

		for z in xrange(info['resolution'][2]):
			for y in xrange(info['resolution'][1]):
				output += {True:' ',False:'\n'}[options.twoD].join(map(lambda x: ('%%%ii'%formatwidth)%(round(x)), distance[i,:,y,z])) + '\n'
			
    
# ------------------------------------------ output result ---------------------------------------  

		file['output'][i].write(output)
	
		if file['name'] != 'STDIN':
			file['output'][i].close()
    
  if file['name'] != 'STDIN':
    file['input'].close()                                                   # close input geom file

    
